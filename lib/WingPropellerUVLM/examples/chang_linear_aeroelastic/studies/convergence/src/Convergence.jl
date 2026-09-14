module Convergence
using Dates, DelimitedFiles, SHA, TOML
module Rigid
include("aerodynamic_backend.jl")
end
const Integrity = Rigid.ChangAerodynamicStudy
module Damping
include("damping_metrics.jl")
end
const FIELDS = (:wing_span,:wing_chord,:prop_radial,:prop_chord,:wake_revolutions,:core,:azimuth_deg)
const MESH = FIELDS[1:4]
const SOURCE = () -> Integrity.source_fingerprint()
const EXAMPLE = normpath(joinpath(@__DIR__,"..","..",".."))

data(x::NamedTuple) = Dict(string(k)=>data(v) for (k,v) in pairs(x))
data(x::AbstractDict) = Dict(string(k)=>data(v) for (k,v) in pairs(x))
data(x::Union{Tuple,AbstractArray}) = [data(v) for v in x]
data(x::Symbol) = string(x)
data(::Nothing) = "nothing"
data(x) = x
named(d::AbstractDict) = (; (Symbol(k)=>v for (k,v) in d)...)
function toml_text(x)
    io=IOBuffer(); TOML.print(io,data(x);sorted=true); return String(take!(io))
end
digest(x) = bytes2hex(sha256(toml_text(x)))
function save_toml(path,x)
    mkpath(dirname(path)); write(path*".tmp",toml_text(x)); mv(path*".tmp",path;force=true)
    return path
end
controls(c) = NamedTuple{FIELDS}(Tuple(getproperty(c,k) for k in FIELDS))
key(c) = Tuple(controls(c))
case(family,level,c) = merge(controls(c),(;family,level,label="$(family)_L$level"))
decreasing(k) = k in (:core,:azimuth_deg)
function validate_controls(c)
    all(k -> hasproperty(c,k),FIELDS) || error("Nominal controls must specify $(FIELDS)")
    all(k -> getproperty(c,k) isa Integer && !(getproperty(c,k) isa Bool) && getproperty(c,k)>0,MESH) ||
        error("Panel counts must be positive integers")
    all(k -> isfinite(getproperty(c,k)) && getproperty(c,k)>0,FIELDS) || error("Controls must be finite and positive")
    n=360/c.azimuth_deg
    n>=4 && isapprox(n,round(n);atol=1e-10,rtol=0) || error("Azimuth step must divide 360 and be <=90 degrees")
end
function build_cases(nominal,sweeps,families; elastic=false)
    validate_controls(nominal)
    !isempty(families) && length(unique(families))==length(families) || error("Families must be nonempty and unique")
    result=NamedTuple[]
    for family in families
        family in FIELDS || error("Unknown family $family")
        elastic && family==:wing_span && error("Wing span also sets the structural mesh; keep it fixed during damping convergence")
        hasproperty(sweeps,family) || error("Missing sweep $family")
        levels=getproperty(sweeps,family)
        length(levels)>=3 || error("$family needs at least three ordered levels")
        all(x -> isfinite(x) && x>0,levels) || error("Invalid $family levels")
        all(x -> decreasing(family) ? x<0 : x>0,diff(levels)) || error("$family levels must be strictly ordered")
        elastic && first(levels)!=getproperty(nominal,family) && error("$family must start at the aerodynamic selection")
        for (level,value) in enumerate(levels)
            c=case(family,level,merge(nominal,NamedTuple{(family,)}((value,))))
            validate_controls(c); push!(result,c)
        end
    end
    return result
end
function core_factors(mode,value)
    mode in (:fixed,:chord,:segment) || error("core_mode must be :fixed, :chord or :segment")
    return (;radius=mode==:fixed ? value : nothing,chord=mode==:chord ? value : 0.,segment=mode==:segment ? value : 0.)
end
function aerodynamic_options(p,mode,c,s)
    f=core_factors(mode,c.core)
    options=Integrity.ChangCoupledAerodynamicOptions(
        flow_speed_mps=p.speed_mps,air_density_kgpm3=p.density_kgpm3,
        angle_of_attack_deg=p.alpha_deg,sideslip_deg=p.beta_deg,
        wing_span_m=p.wing_span_m,wing_root_chord_m=p.wing_root_chord_m,wing_tip_chord_m=p.wing_tip_chord_m,
        wing_symmetric=p.wing_symmetric,wing_span_panels=c.wing_span,wing_chord_panels=c.wing_chord,
        elastic_axis_fraction=p.elastic_axis_fraction,propeller_radius_m=p.propeller_radius_m,
        propeller_chord_m=p.propeller_chord_m,propeller_blades=p.blades,
        propeller_radial_panels=c.prop_radial,propeller_chord_panels=c.prop_chord,
        propeller_attachment_eta=p.attachment_eta,pylon_length_m=p.pylon_length_m,
        reference_rpm=p.rpm,reference_speed_mps=p.speed_mps,collective_pitch_offset_deg=p.collective_offset_deg,
        interaction_on=p.interaction,wake_relaxation=p.wake_shedding_fraction,
        azimuth_step_deg=c.azimuth_deg,retained_wake_revolutions=c.wake_revolutions,
        simulated_revolutions=s.simulated_revolutions,averaged_revolutions=s.averaged_revolutions,
        core_radius_m=f.radius,finite_core_segment_factor=f.segment,finite_core_chord_factor=f.chord)
    Integrity.validate_options(options)
    return options
end
function validate_aerodynamic(s)
    all(x -> isfinite(x) && x>=0,(values(s.absolute_tolerances)...,s.relative_tolerance,values(s.periodic_tolerances)...)) || error("Invalid CT/CQ tolerances")
    cases=build_cases(s.nominal,s.sweeps,s.families)
    for c in cases; aerodynamic_options(s.physical,s.core_mode,c,s); end
    return cases
end
function read_aerodynamic(directory,options,s)
    Integrity.output_matches_options(directory,options;warn_on_mismatch=false) || error("Aerodynamic metadata/options/source/history mismatch: $directory")
    history=joinpath(directory,"coupled_aerodynamic_history.csv")
    raw,header=readdlm(history,',',header=true)
    names=String.(vec(header))
    column(name)=Float64.(raw[:,something(findfirst(==(name),names))])
    steps=Integrity.validate_options(options)
    n=steps*options.simulated_revolutions
    size(raw,1)==n || error("Incomplete aerodynamic history")
    dt=options.azimuth_step_deg/(6*options.reference_rpm*options.flow_speed_mps/options.reference_speed_mps)
    all(isapprox.(column("time_s"),(1:n).*dt;rtol=1e-10,atol=1e-12)) || error("Invalid aerodynamic time grid")
    expected=mod.((1:n).*options.azimuth_step_deg,360)
    all(abs.(mod.(column("azimuth_deg").-expected.+180,360).-180).<1e-8) || error("Invalid aerodynamic phase grid")
    metrics=map(name -> Integrity.periodic_metrics(column(name),steps,options.averaged_revolutions),("propeller_CT","propeller_CQ"))
    return (;status="completed",reason="CT/CQ history verified",values=map(m->m.mean,metrics),
        phases=map(m->m.phase,metrics),periodic=map(m->m.periodic_rms,metrics),
        valid=all(i->metrics[i].periodic_rms<=s.periodic_tolerances[i],1:2),history,
        log=joinpath(directory,"run.log"))
end
function aerodynamic_case(c,s,directory)
    options=aerodynamic_options(s.physical,s.core_mode,c,s); mkpath(directory)
    reused=s.reuse_existing && Integrity.output_matches_options(directory,options;warn_on_mismatch=false)
    if !reused
        open(joinpath(directory,"run.log"),"w") do io
            redirect_stdout(io) do
                redirect_stderr(io) do
                    result=Rigid.simulate_coupled_aerodynamics(options)
                    Rigid.write_results(result,directory;plot_results=false)
                end
            end
        end
    end
    return merge(read_aerodynamic(directory,options,s),(;reused))
end
function limits(reference,s,elastic)
    elastic && return (s.damping.lambda_tolerance_per_s,s.damping.lambda_tolerance_per_s,
        s.damping.frequency_tolerance_hz,s.damping.frequency_tolerance_hz)
    return ntuple(i->max(s.absolute_tolerances[i],s.relative_tolerance*abs(reference.values[i])),2)
end
function errors(a,b,elastic)
    return [elastic ? abs(a.values[i]-b.values[i]) :
        max(abs(a.values[i]-b.values[i]),Integrity.phase_rms_error(a.phases[i],b.phases[i])) for i in eachindex(a.values)]
end
agrees(a,b,s,elastic=false) = a.valid && b.valid && all(errors(a,b,elastic).<=collect(limits(b,s,elastic)))
function annotate(results,cases,s;elastic=false)
    rows=NamedTuple[]
    for family in unique(c.family for c in cases)
        group=sort(filter(r->r.case.family==family,results);by=r->r.case.level)
        isempty(group) && continue
        ref=last(group); expected=count(c->c.family==family,cases)
        complete=length(group)==expected && [r.case.level for r in group]==collect(1:expected) && all(r->r.status=="completed",group)
        stable=complete && length(group)>=2 && agrees(group[end-1],ref,s,elastic)
        for (i,r) in enumerate(group)
            err=r.valid && ref.valid ? errors(r,ref,elastic) : fill(NaN,length(r.values))
            previous=i>1 && r.valid && group[i-1].valid ? errors(r,group[i-1],elastic) : fill(NaN,length(r.values))
            push!(rows,merge(r,(;complete,reference=ref.case.label,reference_stable=stable,
                accepted=stable && all(x->agrees(x,ref,s,elastic),group[i:end]),errors=err,previous_errors=previous,limits=limits(ref,s,elastic))))
        end
    end
    return rows
end
metric_names(elastic) = elastic ? ("pitch_lambda_per_s","yaw_lambda_per_s","pitch_frequency_hz","yaw_frequency_hz") : ("CT","CQ")
function write_tables(directory,rows,cases,p,mode;elastic=false)
    mkpath(directory)
    open(joinpath(directory,"case_matrix.csv"),"w") do io
        writedlm(io,permutedims(["family","level","label",string.(FIELDS)...,"core_mode","dt_s","wake_rows"]),',')
        for c in cases
            writedlm(io,permutedims(Any[c.family,c.level,c.label,Tuple(controls(c))...,mode,
                c.azimuth_deg/(6p.rpm),ceil(Int,c.wake_revolutions*360/c.azimuth_deg)]),',')
        end
    end
    names=metric_names(elastic)
    open(joinpath(directory,"convergence.csv"),"w") do io
        headers=["family","level","label","status","reason","reused","valid","family_complete","reference_stable","accepted","reference",
            names..., ["error_$n" for n in names]..., ["previous_change_$n" for n in names]...,
            ["limit_$n" for n in names]..., (elastic ? String[] : ["periodic_CT","periodic_CQ"])...,"elapsed_s","history_path","log_path"]
        writedlm(io,permutedims(headers),',')
        for r in rows
            writedlm(io,permutedims(Any[r.case.family,r.case.level,r.case.label,r.status,r.reason,r.reused,r.valid,r.complete,
                r.reference_stable,r.accepted,r.reference,r.values...,r.errors...,r.previous_errors...,r.limits...,r.periodic...,
                r.elapsed_s,r.history,r.log]),',')
        end
    end
end
function plot_results(directory,rows;elastic=false)
    panels=Any[]; names=metric_names(elastic)
    for family in unique(r.case.family for r in rows), (i,name) in enumerate(names)
        group=sort(filter(r->r.case.family==family,rows);by=r->r.case.level)
        p=Plots.plot([r.case.level for r in group],[r.values[i] for r in group];marker=:circle,label=name,
            title=string(family),xlabel="Level (case_matrix.csv)",ylabel=name,left_margin=8Plots.mm,bottom_margin=6Plots.mm,
            xticks=[r.case.level for r in group])
        push!(panels,p)
    end
    isempty(panels) && return
    Plots.savefig(Plots.plot(panels...;layout=(cld(length(panels),2),2),size=(1200,320cld(length(panels),2))),joinpath(directory,"convergence.png"))
end
function execute_cases(cases,s,p,mode;elastic=false,runner)
    directory=abspath(s.output_directory); mkpath(directory)
    fingerprint=SOURCE(); save_toml(joinpath(directory,"study.toml"),(;source_sha256=fingerprint,settings=s))
    write_tables(directory,NamedTuple[],cases,p,mode;elastic)
    println("$(elastic ? "Aeroelastic" : "Aerodynamic"): $(length(cases)) entries, $(length(unique(key.(cases)))) unique cases")
    println("Case matrix: $(joinpath(directory,"case_matrix.csv"))")
    s.dry_run && return NamedTuple[]
    results=NamedTuple[]; cache=Dict{Tuple,NamedTuple}()
    for (index,c) in enumerate(cases)
        SOURCE()==fingerprint || error("Solver source changed during study")
        println("[$index/$(length(cases))] $(c.label): $(controls(c))")
        start=time(); folder=joinpath(directory,c.label)
        result=if haskey(cache,key(c))
            merge(cache[key(c)],(;reused=true))
        else
            r=try
                runner(c,folder)
            catch exception
                exception isa InterruptException && rethrow()
                @warn "Case failed" label=c.label exception
                n=elastic ? 4 : 2
                (;status="failed",reason=sprint(showerror,exception),values=ntuple(_->NaN,n),
                    phases=elastic ? () : (Float64[],Float64[]),periodic=elastic ? () : (NaN,NaN),
                    valid=false,history="",log=joinpath(folder,"run.log"),reused=false)
            end
            cache[key(c)]=r; r
        end
        SOURCE()==fingerprint || error("Solver source changed during a case")
        push!(results,merge(result,(;case=c,elapsed_s=time()-start)))
        write_tables(directory,annotate(results,cases,s;elastic),cases,p,mode;elastic)
        println("  $(result.status); valid=$(result.valid); values=$(result.values)")
    end
    return results
end

include("aerodynamic_selection.jl")
include("aeroelastic_adapter.jl")
end
