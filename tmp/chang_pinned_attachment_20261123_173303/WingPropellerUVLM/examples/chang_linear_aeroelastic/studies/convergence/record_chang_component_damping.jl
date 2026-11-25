# Apply the recorded moving-block settings to wing sweeps, using the saved
# wing-tip histories. Blade, wake and azimuth sweeps use pylon pitch and yaw.
import WingPropellerUVLMStudies as Studies
include(joinpath(@__DIR__,"plot_chang_convergence_sensitivity.jl"))
const Convergence = Studies.Convergence

component_target(family) = family in (:wing_span,:wing_chord) ? :wing : :propeller
decode_optional(value) = value == "nothing" ? nothing : value

function component_fit_settings(settings,target)
    d = Convergence.named(settings["damping"])
    r = Convergence.named(settings["response"])
    damping = merge(d,(; target,fit_start_s=decode_optional(d.fit_start_s),
        block_size_mode=Symbol(d.block_size_mode),frequency_band_hz=Tuple(d.frequency_band_hz)))
    response = merge(r,(; impulse_start_s=decode_optional(r.impulse_start_s)))
    integration = Convergence.named(settings["integration"])
    return (; damping,response,integration,make_plots=false)
end

function component_history(directory,row,matrix)
    planned = matrix[row["label"]]
    label = basename(dirname(row["history_path"]))
    source = matrix[label]
    all(k -> planned[k] == source[k],Saved.CONTROLS) ||
        error("Reused history has different controls for $(row["label"])")
    folder = joinpath(directory,label)
    occursin(r"(?m)^validation_passed = true\r?$",read(joinpath(folder,"response_summary.txt"),String)) ||
        error("Incomplete or invalid solver response in $folder")
    return joinpath(folder,"response_history.csv")
end

function record_component_cases(directory,settings)
    planned = Saved.read_table(joinpath(directory,"case_matrix.csv"))
    matrix = Dict(r["label"] => r for r in planned)
    raw = Saved.read_table(joinpath(directory,"convergence.csv"))
    length(raw) == length(planned) || error("The run must have a record for each planned case")
    cases = [(; family=Symbol(r["family"]),level=Int(Saved.number(r["level"])),label=string(r["label"]))
        for r in planned]
    archived = Dict(r["label"] => r for r in raw)
    modes = Dict(target => Studies.MovingBlockStudy(component_fit_settings(settings,target))
        for target in (:wing,:propeller))
    operating = Convergence.named(settings["physical"])
    savedTarget = Symbol(get(settings["damping"],"target","propeller"))
    sameEstimator = get(settings,"damping_method_sha256","") == modes[savedTarget].estimator_sha256
    cache = Dict{Tuple{String,Symbol},NamedTuple}()
    records = mktempdir() do temp
        map(cases) do case
            row = archived[case.label]
            history = component_history(directory,row,matrix)
            target = component_target(case.family)
            metrics = target == :wing ? WING_METRICS : Saved.METRICS
            channels = Convergence.Damping.response_channels(target)
            result = if target == savedTarget && sameEstimator
                fit = TOML.parsefile(joinpath(dirname(history),"damping_fit.toml"))
                (; values=Tuple(Saved.number(row[k]) for k in metrics),valid=Saved.truth(row["valid"]),
                    fit_r_squared=Tuple(fit[string(name)]["fit_r_squared"] for name in channels.names),
                    source="saved fit")
            else
                get!(cache,(history,target)) do
                    input = joinpath(temp,"response_history.csv")
                    cp(history,input;force=true)
                    fit = Convergence.damping_result(input,modes[target],operating)
                    a,b = (getproperty(fit,name) for name in channels.names)
                    (; values=(a.moving_block_lambda_per_s,b.moving_block_lambda_per_s,a.frequency_hz,b.frequency_hz),
                        valid=fit.valid,fit_r_squared=(a.fit_r_squared,b.fit_r_squared),source="refit saved history")
                end
            end
            status = result.valid ? "completed" : "indeterminate"
            (; case,target,history,status,result...)
        end
    end
    println("Refitted $(length(cache)) distinct saved histories")
    return Convergence.annotate(records,cases,modes[:wing];elastic=true),matrix,modes[:wing].estimator_sha256
end

function record_damping_by_component(directory)
    document = TOML.parsefile(joinpath(directory,"study.toml"))
    settings = document["settings"]
    settings["damping_method"] == "moving_block_2" || error("Use an archived moving-block study")
    rows,matrix,fingerprint = record_component_cases(directory,settings)
    output = mkpath(joinpath(directory,"plots","component"))
    open(joinpath(output,"case_damping.csv"),"w") do io
        println(io,"family,level,label,target,mesh_or_step,first_response,second_response,first_lambda_per_s,second_lambda_per_s,first_frequency_hz,second_frequency_hz,first_fit_r_squared,second_fit_r_squared,valid,reference_stable,accepted,reference,first_reference_difference_per_s,second_reference_difference_per_s,fit_source,history_path")
        for r in rows
            family = string(r.case.family)
            signals = r.target == :wing ? WING_SIGNALS : PROPELLER_SIGNALS
            fields = (family,r.case.level,r.case.label,r.target,matrix[r.case.label][family],signals...,
                r.values...,r.fit_r_squared...,r.valid,r.reference_stable,r.accepted,r.reference,
                r.errors[1],r.errors[2],r.source,relpath(r.history,directory))
            writedlm(io,permutedims(collect(fields)),',')
        end
    end
    Convergence.save_toml(joinpath(output,"damping_settings.toml"),
        (; history_source_sha256=document["source_sha256"],damping_method_sha256=fingerprint,
            response=settings["response"],damping=settings["damping"],
            wing_families=["wing_span","wing_chord"],propeller_families=["prop_radial","prop_chord","wake_revolutions","azimuth_deg"]))
    groups = map(AXIS_LABELS) do (family,xlabel)
        selected = filter(r -> string(r.case.family) == family,rows)
        target = component_target(Symbol(family))
        signals = target == :wing ? WING_SIGNALS : PROPELLER_SIGNALS
        points = map(eachindex(selected)) do i
            r = selected[i]
            x = Saved.number(matrix[r.case.label][family])
            previous = i == 1 ? ntuple(_ -> NaN,4) : abs.(r.values .- selected[i-1].values)
            successive_percent = ntuple(j -> Saved.relative_error_percent(previous[j],r.values[j]),4)
            successive_valid = r.valid && i > 1 && selected[i-1].valid
            (; x,r.values,r.valid,previous,successive_percent,successive_valid)
        end
        return (; family,xlabel,signals,rows=points)
    end
    return save_sensitivity_plots(collect(groups),output,settings["physical"]["speed_mps"])
end

if abspath(PROGRAM_FILE) == @__FILE__
    println(record_damping_by_component(abspath(only(ARGS))))
end
