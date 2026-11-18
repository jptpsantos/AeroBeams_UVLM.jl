# Plot archived damping results without loading a study or rerunning simulations.
# julia --project=lib/WingPropellerUVLM/studies path/to/this_script.jl [results_directory]
# Add --current to plot the latest root run instead of the two archived cases.
# Add --mesh-case1 to combine Mesh_case1 with the latest root run (or Mesh_case2).
# Add --reference-style for compact serif figures: pitch above yaw, mesh directions across.
module ChangConvergenceAeroelasticPlots
using DelimitedFiles,Plots,Printf,TOML

const DEFAULT_DIRECTORY = normpath(joinpath(@__DIR__,"..","..","output","convergence_aeroelastic_explicit"))
const CASE_NAMES = ("First case","Second case")
const FAMILIES = ("prop_radial","prop_chord")
const CONTROLS = ("wing_span","wing_chord","prop_radial","prop_chord",
    "wake_revolutions","core","azimuth_deg")
const METRICS = ("pitch_lambda_per_s","yaw_lambda_per_s","pitch_frequency_hz","yaw_frequency_hz")
const COLORS = (:steelblue,:chocolate3)
const MARKERS = (:circle,:diamond)

number(x) = x isa Number ? Float64(x) : parse(Float64,string(x))
truth(x) = lowercase(string(x)) == "true"

function read_table(path)
    raw,header = readdlm(path,','; header=true)
    columns = Dict(string(name) => i for (i,name) in enumerate(vec(header)))
    return [Dict(name => raw[i,column] for (name,column) in columns) for i in axes(raw,1)]
end

# Retain fingerprints and all settings that affect the response or estimator.
function comparison_settings(document)
    comparable = deepcopy(document)
    for key in ("families","sweeps","output_directory","dry_run","make_plots")
        pop!(comparable["settings"],key,nothing)
    end
    return comparable
end

same_numbers(a,b) = all(isapprox.(a,b; rtol=1e-10,atol=1e-12))

# A relative error is undefined at zero growth/decay rate.
relative_error_percent(error,value) =
    !isfinite(error) || !isfinite(value) || iszero(value) ? NaN : 100error / abs(value)
successive_errors_percent(row) = ntuple(i -> relative_error_percent(row.previous[i],row.values[i]),2)
relative_limits_percent(row) = ntuple(i -> relative_error_percent(row.limits[i],row.values[i]),2)
reference_errors_percent(row,reference) =
    ntuple(i -> relative_error_percent(row.errors[i],reference.values[i]),2)

"""Combine the two archived blade sweeps using full controls as the case identity.

`cases` contains each physical/numerical configuration once. `rows` contains its
membership in each sweep, so the common baseline can be used in both panels.
Saved errors and level numbers are discarded and recomputed on the merged grid.
"""
function load_results(directory=DEFAULT_DIRECTORY; case_names=CASE_NAMES,
    families=FAMILIES,allow_invalid=false)
    isempty(case_names) && error("At least one archived case is required")
    cases = NamedTuple[]
    memberships = Dict{Tuple{String,Float64},Int}()
    identities = Dict{Tuple,Int}()
    comparison = nothing
    settings = nothing
    input_count = 0
    for caseName in case_names
        archive = joinpath(directory,caseName)
        document = TOML.parsefile(joinpath(archive,"study.toml"))
        comparable = comparison_settings(document)
        if isnothing(comparison)
            comparison = comparable
            settings = document["settings"]
        else
            isequal(comparison,comparable) || error("Incompatible saved study settings or fingerprints: $caseName")
        end
        matrix = Dict(
            (string(c["family"]),Int(number(c["level"]))) => c
            for c in read_table(joinpath(archive,"case_matrix.csv"))
        )
        saved = filter(r -> string(r["family"]) in families,
            read_table(joinpath(archive,"convergence.csv")))
        isempty(saved) && error("No blade convergence results in $archive")
        for r in saved
            input_count += 1
            family = string(r["family"])
            c = matrix[(family,Int(number(r["level"])))]
            label = string(r["label"])
            string(c["label"]) == label || error("Case label mismatch: $caseName/$label")
            valid = string(r["status"]) == "completed" && truth(r["valid"])
            valid || allow_invalid ||
                error("Cannot compare an incomplete or invalid case: $caseName/$label")
            controls = Tuple(number(c[k]) for k in CONTROLS)
            all(x -> isfinite(x) && x > 0,controls) || error("Invalid controls: $caseName/$label")
            all(isinteger,controls[1:4]) || error("Noninteger panel counts: $caseName/$label")
            for (k,value) in zip(CONTROLS,controls)
                k == family && continue
                value == number(settings["nominal"][k]) ||
                    error("Non-swept control $k differs from nominal: $caseName/$label")
            end
            string(c["core_mode"]) == settings["core_mode"] || error("Core law mismatch: $caseName/$label")
            isapprox(number(c["dt_s"]),controls[7] / (6settings["physical"]["rpm"]); rtol=1e-10) ||
                error("Time step mismatch: $caseName/$label")
            values = Tuple(number(r[k]) for k in METRICS)
            limits = Tuple(number(r["limit_$k"]) for k in METRICS)
            all(isfinite,values) && all(x -> isfinite(x) && x > 0,limits) ||
                error("Invalid damping results or tolerances: $caseName/$label")
            all(>(0),values[3:4]) || error("Invalid fitted frequency: $caseName/$label")
            d = settings["damping"]
            expectedLimits = (d["lambda_tolerance_per_s"],d["lambda_tolerance_per_s"],
                d["frequency_tolerance_hz"],d["frequency_tolerance_hz"])
            same_numbers(limits,expectedLimits) || error("Tolerance mismatch: $caseName/$label")
            identity = (controls...,string(c["core_mode"]))
            source = "$caseName/$label"
            index = get(identities,identity,0)
            if index == 0
                index = length(cases) + 1
                identities[identity] = index
                push!(
                    cases,
                    (; id=@sprintf("case_%03d",index),controls,values,limits,
                        valid,status=string(r["status"]),reason=string(r["reason"]),
                        sources=String[source],families=String[family])
                )
            else
                canonical = cases[index]
                same_numbers(values,canonical.values) && same_numbers(limits,canonical.limits) ||
                    error("Conflicting duplicate results: $source and $(first(canonical.sources))")
                canonical.valid == valid || error("Conflicting duplicate validity: $source")
                push!(canonical.sources,source)
                family in canonical.families || push!(canonical.families,family)
            end
            member = (family,number(c[family]))
            get(memberships,member,index) == index ||
                error("Different configurations at the same $family panel count")
            memberships[member] = index
        end
    end
    rows = NamedTuple[]
    for family in families
        levels = sort([x for (f,x) in keys(memberships) if f == family])
        length(levels) >= 3 || error("$family needs at least three distinct panel counts")
        reference = cases[memberships[(family,last(levels))]]
        for (i,x) in enumerate(levels)
            c = cases[memberships[(family,x)]]
            previous = i == 1 ? ntuple(_ -> NaN,4) :
                       abs.(c.values .- cases[memberships[(family,levels[i-1])]].values)
            errors = abs.(c.values .- reference.values)
            push!(
                rows,
                (; family,x,case_id=c.id,values=c.values,limits=c.limits,
                    previous,errors,reference_panels=last(levels))
            )
        end
    end
    return (; rows,cases,settings,input_count)
end

function panel_options(family,group,nominal)
    radial = family == "prop_radial"
    ticks = [r.x for r in group]
    padding = 0.04(last(ticks) - first(ticks))
    fixedPanels = nominal[radial ? "prop_chord" : "prop_radial"]
    return (;
        title=radial ? "Blade radial panels (chord = $fixedPanels)" : "Blade chord panels (radial = $fixedPanels)",
        xlabel=radial ? "Radial panels per blade" : "Chordwise panels per blade",
        xticks=(ticks,[@sprintf("%g",x) for x in ticks]),
        xlims=(first(ticks) - padding,last(ticks) + padding),
        titlefontsize=11,guidefontsize=10,tickfontsize=9,
        legendfontsize=9,gridalpha=0.16,framestyle=:box,
        left_margin=8Plots.mm,bottom_margin=7Plots.mm,
        top_margin=4Plots.mm,right_margin=4Plots.mm)
end

function damping_panels(rows,nominal)
    panels = Any[]
    lo,hi = extrema([r.values[i] for r in rows for i = 1:2])
    padding = max(0.08(hi - lo),0.01max(abs(lo),abs(hi)))
    errorValues = [max(value,eps(Float64)) for r in rows for value in successive_errors_percent(r) if isfinite(value)]
    toleranceValues = [value for r in rows for value in relative_limits_percent(r) if isfinite(value)]
    isempty(errorValues) && error("No finite relative damping errors to plot")
    errorLo = min(minimum(errorValues),minimum(toleranceValues)) / 1.5
    errorHi = max(maximum(errorValues),maximum(toleranceValues)) * 1.5
    errorTicks = filter(x -> errorLo <= x <= errorHi,
        [m * 10.0^e for e = floor(Int,log10(errorLo)):ceil(Int,log10(errorHi)) for m in (1,2,5)])
    for family in FAMILIES
        group = filter(r -> r.family == family,rows)
        p = plot(; panel_options(family,group,nominal)...,ylabel="Damping rate, λ (s⁻¹)",
            ylims=(lo - padding,hi + padding),legend=:bottomright)
        for (i,name) in enumerate(("Pitch","Yaw"))
            plot!(p,[r.x for r in group],[r.values[i] for r in group];
                label=name,color=COLORS[i],marker=MARKERS[i],
                markersize=i == 1 ? 5 : 3.5,markerstrokecolor=COLORS[i],
                linestyle=i == 1 ? :solid : :dash,linewidth=2)
        end
        push!(panels,p)
    end
    for family in FAMILIES
        group = filter(r -> r.family == family,rows)
        refined = group[2:end]
        p = plot(; panel_options(family,group,nominal)...,
            ylabel="Successive relative damping error (%)",
            yscale=:log10,ylims=(errorLo,errorHi),legend=:topright,
            yticks=(errorTicks,[@sprintf("%g",x) for x in errorTicks]))
        for (i,name) in enumerate(("Pitch","Yaw"))
            plot!(p,[r.x for r in refined],[max(successive_errors_percent(r)[i],eps(Float64)) for r in refined];
                label=name,color=COLORS[i],marker=MARKERS[i],
                markersize=i == 1 ? 5 : 3.5,markerstrokecolor=COLORS[i],
                linestyle=i == 1 ? :solid : :dash,linewidth=2)
        end
        for (i,name) in enumerate(("Pitch","Yaw"))
            plot!(p,[r.x for r in refined],[relative_limits_percent(r)[i] for r in refined];
                label="$name tolerance (%)",color=COLORS[i],linestyle=:dot,
                linewidth=1.3,alpha=0.7)
        end
        push!(panels,p)
    end
    return panels
end

"""Show successive and finest-reference damping differences on linear percent axes."""
function damping_difference_panels(rows,nominal)
    panels = Any[]
    for comparison in (:successive,:reference)
        curves = map(FAMILIES) do family
            group = sort(filter(r -> r.family == family,rows); by=r -> r.x)
            reference = last(group)
            points = comparison == :successive ? group[2:end] : group
            errors = [
                comparison == :successive ? successive_errors_percent(r) :
                reference_errors_percent(r,reference) for r in points
            ]
            limits = [relative_limits_percent(comparison == :successive ? r : reference)
                      for r in points]
            (; family,group,reference,points,errors,limits)
        end
        finiteValues = [
            value for curve in curves for pair in (curve.errors...,curve.limits...)
            for value in pair if isfinite(value)
        ]
        ymax = isempty(finiteValues) ? 1.0 : max(1.12maximum(finiteValues),1e-3)
        for curve in curves
            (; family,group,reference,points,errors,limits) = curve
            options = panel_options(family,group,nominal)
            subtitle =
                comparison == :successive ? "Successive cases" :
                @sprintf("Finest reference: %g panels",reference.x)
            options = merge(options,(; title=options.title * "\n" * subtitle))
            ylabel =
                comparison == :successive ? "Successive damping difference (%)" :
                "Difference from finest damping (%)"
            p = plot(; options...,ylabel,ylims=(0.0,ymax),legend=:topright)
            for (i,name) in enumerate(("Pitch","Yaw"))
                plot!(p,[r.x for r in points],[pair[i] for pair in errors];
                    label=name,color=COLORS[i],marker=MARKERS[i],
                    markersize=i == 1 ? 5 : 3.5,markerstrokecolor=COLORS[i],
                    linestyle=i == 1 ? :solid : :dash,linewidth=2)
                plot!(p,[r.x for r in points],[pair[i] for pair in limits];
                    label="$name tolerance (%)",color=COLORS[i],linestyle=:dot,
                    linewidth=1.3,alpha=0.7)
            end
            push!(panels,p)
        end
    end
    return panels
end

function write_tables(out,results)
    open(joinpath(out,"unique_cases.csv"),"w") do io
        writedlm(io,permutedims(["case_id",CONTROLS...,METRICS...,"families","sources"]),',')
        for c in results.cases
            writedlm(io,
                permutedims(Any[c.id,c.controls...,c.values...,
                    join(c.families,"; "),join(c.sources,"; ")]),',')
        end
    end
    open(joinpath(out,"damping_convergence.csv"),"w") do io
        headers = ["family","panels","case_id","reference_panels",METRICS...,
            ("previous_change_$k" for k in METRICS)...,
            ("error_$k" for k in METRICS)...,
            "previous_relative_pitch_percent","previous_relative_yaw_percent",
            "reference_relative_pitch_percent","reference_relative_yaw_percent",
            "limit_pitch_lambda_per_s","limit_yaw_lambda_per_s"]
        writedlm(io,permutedims(headers),',')
        for r in results.rows
            group = filter(a -> a.family == r.family,results.rows)
            relativePrevious = successive_errors_percent(r)
            reference = last(group)
            relativeReference = reference_errors_percent(r,reference)
            writedlm(io,
                permutedims(
                    Any[r.family,r.x,r.case_id,r.reference_panels,
                        r.values...,r.previous...,r.errors...,relativePrevious...,
                        relativeReference...,r.limits[1:2]...]
                ),',')
        end
    end
end

function write_notes(out,results,case_names)
    s = results.settings
    open(joinpath(out,"plot_notes.md"),"w") do io
        println(io,"# Combined Chang blade-mesh damping convergence\n")
        println(
            io,
            "Sources: $(join(case_names," + ")). Only archived convergence.csv, case_matrix.csv and study.toml are read; no simulations run."
        )
        println(
            io,
            "Source fingerprints and response, structural, physical, integration, coupling and damping settings match."
        )
        println(
            io,
            "Operating point: U = $(s["physical"]["speed_mps"]) m/s; alpha = $(s["physical"]["alpha_deg"]) deg; RPM = $(s["physical"]["rpm"])."
        )
        println(
            io,
            "Fixed wing mesh: $(s["nominal"]["wing_span"]) span x $(s["nominal"]["wing_chord"]) chord; azimuth increment = $(s["nominal"]["azimuth_deg"]) deg; wake = $(s["nominal"]["wake_revolutions"]) rev."
        )
        println(
            io,
            "Core law = $(s["core_mode"]); core value = $(s["nominal"]["core"]). Wing/propeller interaction = $(s["physical"]["interaction"]).\n"
        )
        println(
            io,
            "$(results.input_count) input blade-family rows become $(length(results.rows)) distinct sweep points and $(length(results.cases)) unique configurations."
        )
        println(
            io,
            "unique_cases.csv stores each configuration once with all source labels. Configurations shared by both curves use one case_id."
        )
        println(
            io,
            "Duplicate identity uses all seven numerical controls and the core law; conflicting duplicate results cause an error. Saved level labels are not identities.\n"
        )
        println(
            io,
            "The plotted damping metric is the saved moving-block lambda (1/s); negative lambda means decay. All selected fits are completed and valid."
        )
        println(
            io,
            "Plotted successive error (%) = 100 * abs(lambda(N_i) - lambda(N_(i-1))) / abs(lambda(N_i)) after sorting the combined panel counts. The first point has no previous error."
        )
        println(
            io,
            "Dotted pitch/yaw tolerance curves (%) = 100 * saved absolute lambda tolerance / abs(current lambda). The saved acceptance criterion is unchanged."
        )
        println(
            io,
            "Zero errors are displayed at machine epsilon on the log axis. Relative errors at zero lambda are undefined and saved as NaN."
        )
        println(
            io,
            "The CSV also contains errors against each merged finest-panel reference, recalculated rather than copied from the archives."
        )
        println(
            io,
            "Successive relative error (%) = 100 * successive absolute error / abs(current lambda); reference relative error uses abs(finest lambda).\n"
        )
        println(
            io,
            "damping_percentage_differences.png/.svg shows both comparisons on linear percentage axes: successive cases in the top row and each case versus its family's finest reference in the bottom row. Radial/chord sweeps occupy the left/right columns."
        )
        println(
            io,
            "The finest reference has exactly zero reference difference, displayed at zero on this linear axis. Dotted tolerance curves use the same denominator as the corresponding difference; no simulation acceptance criteria are changed.\n"
        )
        println(
            io,
            "| Sweep | Combined panels | First count whose remaining successive pitch/yaw changes meet tolerance | First count whose remaining pitch/yaw reference errors meet tolerance |"
        )
        println(io,"|---|---|---|---|")
        for family in FAMILIES
            group = filter(r -> r.family == family,results.rows)
            successive = findfirst(
                i -> i > 1 && all(r -> all(r.previous[1:2] .<= r.limits[1:2]),group[i:end]),
                eachindex(group)
            )
            reference = findfirst(i -> all(r -> all(r.errors[1:2] .<= r.limits[1:2]),group[i:end]),eachindex(group))
            format_count(i) = isnothing(i) ? "none" : @sprintf("%g",group[i].x)
            println(
                io,
                "| $family | $(join([@sprintf("%g",r.x) for r in group],", ")) | $(format_count(successive)) | $(format_count(reference)) |"
            )
        end
        println(
            io,
            "\nThese are separate one-parameter mesh sensitivity sweeps at fixed time step, not a test with both blade directions refined together. The finest mesh is a numerical reference, not an exact solution."
        )
        println(
            io,
            "The root directory's separate azimuth convergence run is excluded. These figures do not establish time-step convergence or a flutter boundary."
        )
    end
end

function plot_saved_convergence(directory=DEFAULT_DIRECTORY; case_names=CASE_NAMES)
    directory = abspath(directory)
    results = load_results(directory; case_names)
    out = joinpath(directory,"plots")
    mkpath(out)
    panels = damping_panels(results.rows,results.settings["nominal"])
    p = results.settings["physical"]
    title = @sprintf("%s | U = %g m/s, alpha = %g deg",join(case_names," + "),p["speed_mps"],p["alpha_deg"])
    figure = plot(panels...; layout=(2,2),size=(1400,900),plot_title=title,
        plot_titlefontsize=14,background_color=:white,dpi=180)
    paths = [joinpath(out,"damping_convergence" * extension) for extension in (".png",".svg")]
    foreach(path -> savefig(figure,path),paths)
    differenceFigure = plot(damping_difference_panels(results.rows,results.settings["nominal"])...;
        layout=(2,2),size=(1400,950),
        plot_title="Damping convergence: percentage differences | " * title,
        plot_titlefontsize=13,background_color=:white,dpi=180)
    differencePaths = [joinpath(out,"damping_percentage_differences" * extension)
                       for extension in (".png",".svg")]
    foreach(path -> savefig(differenceFigure,path),differencePaths)
    append!(paths,differencePaths)
    write_tables(out,results)
    write_notes(out,results,case_names)
    println("Merged $(results.input_count) archived rows into $(length(results.cases)) unique configurations.")
    foreach(println,paths)
    return paths
end

"""Read the current run without merging older archived studies.

Indeterminate fits remain diagnostic data and retain their validity flags.
Differences use adjacent saved levels, including gaps with nonfinite values.
"""
function load_current_results(directory=DEFAULT_DIRECTORY; metrics=METRICS)
    saved = read_table(joinpath(directory,"convergence.csv"))
    matrix = Dict((string(c["family"]),Int(number(c["level"]))) => c
                  for c in read_table(joinpath(directory,"case_matrix.csv")))
    isempty(saved) && error("No current convergence results in $directory")
    rows = NamedTuple[]
    for family in unique(string(r["family"]) for r in saved)
        group = sort(filter(r -> string(r["family"]) == family,saved);
                     by=r -> number(r["level"]))
        reference = last(group)
        referenceValues = Tuple(number(reference[k]) for k in metrics)
        previousValues = ntuple(_ -> NaN,4)
        previousValid = false
        for r in group
            level = Int(number(r["level"]))
            c = matrix[(family,level)]
            string(c["label"]) == string(r["label"]) || error("Case label mismatch")
            values = Tuple(number(r[k]) for k in metrics)
            limits = Tuple(number(r["limit_$k"]) for k in metrics)
            valid = truth(r["valid"]) && string(r["status"]) == "completed"
            previous = abs.(values .- previousValues)
            errors = abs.(values .- referenceValues)
            successivePercent = ntuple(i -> relative_error_percent(previous[i],values[i]),4)
            referencePercent = ntuple(i -> relative_error_percent(errors[i],referenceValues[i]),4)
            push!(rows,(; family,level,label=string(r["label"]),x=number(c[family]),
                values,limits,valid,status=string(r["status"]),reason=string(r["reason"]),
                previous,errors,successive_percent=successivePercent,
                reference_percent=referencePercent,
                successive_valid=valid && previousValid,
                reference_valid=valid && truth(reference["valid"]) && string(reference["status"]) == "completed",
                reference_label=string(reference["label"]),sources=string(r["label"])))
            previousValues = values
            previousValid = valid
        end
    end
    settings = TOML.parsefile(joinpath(directory,"study.toml"))["settings"]
    return (; rows,settings)
end

mesh_case_names(directory) =
    ("Mesh_case1",isfile(joinpath(directory,"convergence.csv")) ? "." : "Mesh_case2")

"""Merge Mesh_case1 with the latest run by complete numerical configuration.

Compatibility and duplicate checks are shared with the archived plotter.
Invalid fits remain flagged; differences are recomputed on the merged grid.
"""
function load_mesh_results(directory=DEFAULT_DIRECTORY; case_names=mesh_case_names(directory),
    families=("wing_span","wing_chord","prop_radial","prop_chord"))
    combined = load_results(directory; case_names,families,allow_invalid=true)
    identities = Dict(c.id => c for c in combined.cases)
    rows = NamedTuple[]
    for family in unique(r.family for r in combined.rows)
        group = filter(r -> r.family == family,combined.rows)
        reference = identities[last(group).case_id]
        for (level,r) in enumerate(group)
            c = identities[r.case_id]
            previousValid = level > 1 && identities[group[level-1].case_id].valid
            push!(rows,(; r.family,level,label=c.id,r.x,r.values,r.limits,c.valid,c.status,c.reason,
                r.previous,r.errors,
                successive_percent=ntuple(i -> relative_error_percent(r.previous[i],r.values[i]),4),
                reference_percent=ntuple(i -> relative_error_percent(r.errors[i],reference.values[i]),4),
                successive_valid=c.valid && previousValid,
                reference_valid=c.valid && reference.valid,
                reference_label=reference.id,sources=join(c.sources,"; ")))
        end
    end
    return (; rows,combined.settings,input_count=combined.input_count,cases=combined.cases,combined)
end

"""Compact two-by-two figures matching the supplied black/red reference style.

Lines connect computed values directly. Fit validity remains in the data tables.
"""
function plot_reference_style(results,out)
    paths = String[]
    meshes = (("wing",("wing_span","wing_chord"),("Spanwise panels","Chordwise panels")),
              ("propeller",("prop_radial","prop_chord"),("Radial panels","Chordwise panels")))
    for (name,families,xlabels) in meshes
        all(family -> any(r -> r.family == family,results.rows),families) || continue
        for comparison in (:successive,:reference)
            panels = Any[]
            for (i,metric) in enumerate(("Pitch","Yaw")),(family,xlabel) in zip(families,xlabels)
                group = filter(r -> r.family == family,results.rows)
                points = comparison == :successive ? group[2:end] : group
                differences = [comparison == :successive ? r.successive_percent[i] : r.reference_percent[i]
                               for r in points]
                finiteValues = filter(isfinite,differences)
                ymax = isempty(finiteValues) ? 1.0 : max(1.08 * maximum(finiteValues),1e-3)
                xmax = maximum(r.x for r in group)
                # Regular ticks keep the smallest mesh levels from crowding the axis.
                step = xmax <= 15 ? 5 : 10
                ticks = collect(0:step:xmax)
                p = plot([r.x for r in points],differences;
                    color=i == 1 ? :black : "#EF4444",linestyle=i == 1 ? :solid : :dash,
                    linewidth=1.1,label="",legend=false,grid=false,
                    framestyle=:axes,mirror=false,tick_direction=:in,
                    fontfamily="Times",guidefontsize=9,tickfontsize=8,
                    foreground_color_axis=:black,foreground_color_border=:black,
                    foreground_color_text=:black,background_color=:white,
                    ylabel="$metric difference (%)",xlabel=i == 2 ? xlabel : "",
                    xlims=(0.0,xmax),xticks=ticks,ylims=(0.0,ymax),widen=false,
                    left_margin=20Plots.mm,right_margin=7Plots.mm,
                    top_margin=3Plots.mm,bottom_margin=i == 2 ? 10Plots.mm : 1Plots.mm)
                push!(panels,p)
            end
            figure = plot(panels...; layout=(2,2),size=(1000,390),
                background_color=:white,dpi=220)
            for extension in (".svg",".png")
                path = joinpath(out,"$(name)_damping_$(comparison)_percent" * extension)
                savefig(figure,path)
                push!(paths,path)
            end
        end
    end
    isempty(paths) && error("Reference style requires both mesh directions for a wing or propeller")
    return paths
end

"""Plot current-run percent differences; no simulations or fit changes."""
function plot_current_convergence(directory=DEFAULT_DIRECTORY; case_names=nothing,reference_style=false)
    directory = abspath(directory)
    results = isnothing(case_names) ? load_current_results(directory) : load_mesh_results(directory; case_names)
    out = joinpath(directory,"plots",isnothing(case_names) ? "current_run" : "mesh_case1_and_latest")
    mkpath(out)
    labels = Dict("wing_span" => ("Wing span mesh","Spanwise panels"),
        "wing_chord" => ("Wing chord mesh","Chordwise panels"),
        "prop_radial" => ("Propeller radial mesh","Radial panels per blade"),
        "prop_chord" => ("Propeller chord mesh","Chordwise panels per blade"))
    families = unique(r.family for r in results.rows)
    paths = String[]
    if reference_style
        append!(paths,plot_reference_style(results,out))
    end
    comparisons = reference_style ? () : (:successive,:reference)
    for comparison in comparisons
        panels = Any[]
        for family in families
            group = filter(r -> r.family == family,results.rows)
            points = comparison == :successive ? group[2:end] : group
            title,xlabel = get(labels,family,(family,family))
            valid_pair(r) = comparison == :successive ? r.successive_valid : r.reference_valid
            diagnostic = any(r -> !valid_pair(r),points)
            title *= diagnostic ? "\nDiagnostic: unvalidated fits (x)" : ""
            ticks = [r.x for r in group]
            padding = max(0.04 * (maximum(ticks) - minimum(ticks)),0.01 * maximum(abs.(ticks)),1e-8)
            plottedValues = Float64[]
            for r in points,i in 1:2
                difference = comparison == :successive ? r.successive_percent[i] : r.reference_percent[i]
                limit = relative_error_percent(r.limits[i],
                    comparison == :successive ? r.values[i] : last(group).values[i])
                isfinite(difference) && push!(plottedValues,difference)
                isfinite(limit) && push!(plottedValues,limit)
            end
            ymax = isempty(plottedValues) ? 1.0 : max(1.12 * maximum(plottedValues),1e-3)
            p = plot(; title,xlabel,ylabel="Damping difference (%)",ylims=(0.0,ymax),
                xticks=ticks,xlims=(minimum(ticks) - padding,maximum(ticks) + padding),
                legend=:topright,titlefontsize=11,guidefontsize=10,tickfontsize=9,
                legendfontsize=8,gridalpha=0.16,framestyle=:box,
                left_margin=7Plots.mm,bottom_margin=6Plots.mm,top_margin=4Plots.mm)
            for (i,name) in enumerate(("Pitch","Yaw"))
                differences = [comparison == :successive ? r.successive_percent[i] : r.reference_percent[i]
                               for r in points]
                plot!(p,[r.x for r in points],differences; label=name,color=COLORS[i],
                    linewidth=2,linestyle=i == 1 ? :solid : :dash)
                for validated in (true,false)
                    indices = findall(r -> valid_pair(r) == validated,points)
                    isempty(indices) && continue
                    scatter!(p,[points[j].x for j in indices],differences[indices];
                        label="",color=COLORS[i],marker=validated ? MARKERS[i] : :xcross,
                        markersize=i == 1 ? 6 : 4,markerstrokecolor=COLORS[i])
                end
                limits = [relative_error_percent(r.limits[i],
                    comparison == :successive ? r.values[i] : last(group).values[i]) for r in points]
                plot!(p,[r.x for r in points],limits; label="$name tolerance",color=COLORS[i],
                    linestyle=:dot,linewidth=1.2,alpha=0.7)
            end
            push!(panels,p)
        end
        title = comparison == :successive ? "Successive mesh differences" : "Differences from finest mesh in each sweep"
        isnothing(case_names) || (title *= " | Mesh_case1 + latest")
        physical = results.settings["physical"]
        title *= @sprintf(" | U = %g m/s, alpha = %g deg",physical["speed_mps"],physical["alpha_deg"])
        figure = plot(panels...; layout=(cld(length(panels),2),2),
            size=(1400,450 * cld(length(panels),2)),plot_title=title,
            plot_titlefontsize=13,background_color=:white,dpi=180)
        for extension in (".png",".svg")
            path = joinpath(out,"damping_$(comparison)_percent" * extension)
            savefig(figure,path)
            push!(paths,path)
        end
    end
    open(joinpath(out,"percentage_differences.csv"),"w") do io
        headers = ["family","level","label","control_value","status","valid","reason","sources",
            "reference_label","successive_valid","reference_valid",
            ("successive_percent_$k" for k in METRICS)...,
            ("reference_percent_$k" for k in METRICS)...]
        writedlm(io,permutedims(headers),',')
        for r in results.rows
            writedlm(io,permutedims(Any[r.family,r.level,r.label,r.x,r.status,r.valid,r.reason,r.sources,
                r.reference_label,r.successive_valid,r.reference_valid,
                r.successive_percent...,r.reference_percent...]),',')
        end
    end
    isnothing(case_names) || write_tables(out,results.combined)
    open(joinpath(out,"plot_notes.md"),"w") do io
        println(io,"# Damping percentage differences\n")
        println(io,"Source: $directory. Reads saved CSV/TOML files only; no simulations or damping fits run.\n")
        if !isnothing(case_names)
            println(io,"Combined sources: $(join(case_names," + ")) (. means the latest root run). Study settings and source fingerprints match after excluding sweep lists and output controls.")
            println(io,"$(results.input_count) input rows produce $(length(results.rows)) sweep points and $(length(results.cases)) unique configurations. Duplicate configurations are checked for matching metrics, limits and validity.\n")
            println(io,"Panel counts are sorted before recomputing differences. CSV sources identify every contributing case; source level numbers are not configuration identities.")
        end
        println(io,"Successive (%) = 100 * abs(lambda_i - lambda_previous) / abs(lambda_i). The first level is undefined (NaN).")
        println(io,"Finest reference (%) = 100 * abs(lambda_i - lambda_finest) / abs(lambda_finest). The finest reference is the last saved refinement level in each family, not an exact solution.")
        println(io,"Zero denominators and unavailable values give NaN. Dotted tolerances use the saved absolute damping tolerance with the corresponding denominator.\n")
        if reference_style
            println(io,"Compact reference style: serif type, black solid pitch curves above red dashed yaw curves, mesh directions in columns, no grid, legend, tolerance curves or markers. Lines directly connect saved mesh levels; no interpolation is applied.")
            println(io,"Each mesh has separate successive and finest-reference figures. Fit validity is retained in the CSV rather than marked on these compact curves.")
        else
            println(io,"Crosses mark comparisons involving a fit that was not completed and valid.")
        end
        println(io,"These are diagnostic differences, not verified convergence. Saved validity, acceptance flags and source files are unchanged.")
        println(io,"$(count(r -> !r.valid,results.rows)) of $(length(results.rows)) saved mesh points have unvalidated fits. The CSV retains their status and reason.")
        println(io,"The CSV also includes the corresponding frequency differences; the figures show damping rates only.")
    end
    foreach(println,paths)
    return paths
end
end

if abspath(PROGRAM_FILE) == @__FILE__
    current = "--current" in ARGS
    meshCase1 = "--mesh-case1" in ARGS
    referenceStyle = "--reference-style" in ARGS
    directories = filter(arg -> !(arg in ("--current","--mesh-case1","--reference-style")),ARGS)
    directory = isempty(directories) ? ChangConvergenceAeroelasticPlots.DEFAULT_DIRECTORY : only(directories)
    if meshCase1
        ChangConvergenceAeroelasticPlots.plot_current_convergence(directory;
            case_names=ChangConvergenceAeroelasticPlots.mesh_case_names(directory),reference_style=referenceStyle)
    elseif current
        ChangConvergenceAeroelasticPlots.plot_current_convergence(directory; reference_style=referenceStyle)
    else
        ChangConvergenceAeroelasticPlots.plot_saved_convergence(directory)
    end
end
