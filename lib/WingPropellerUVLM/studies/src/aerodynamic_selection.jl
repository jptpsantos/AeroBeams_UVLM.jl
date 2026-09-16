function selected_controls(rows,s)
    Set(s.families)==Set(FIELDS) || error("A selection requires all seven aerodynamic families; partial sweeps only produce inspection results")
    selected=s.nominal
    for family in FIELDS
        group=sort(filter(r->r.case.family==family && r.accepted,rows);by=r->r.case.level)
        isempty(group) && error("No verified CT/CQ plateau for $family; extend/refine that study or its settling duration")
        value=getproperty(first(group).case,family)
        selected=merge(selected,NamedTuple{(family,)}((value,)))
    end
    return controls(selected)
end
function confirmation_cases(selected,s)
    finer=selected
    # Keep the selected core law/value, refining discretization around it.
    for k in FIELDS
        k==:core && continue
        current=getproperty(selected,k)
        candidates=filter(v->decreasing(k) ? v<current : v>current,getproperty(s.sweeps,k))
        isempty(candidates) && error("No finer $k level available for combined verification")
        finer=merge(finer,NamedTuple{(k,)}((first(candidates),)))
    end
    return [case(:combined,1,selected),case(:combined,2,finer)]
end
function selection_payload(s,selected,results)
    return Dict("schema"=>1,"status"=>"verified","source_sha256"=>SOURCE(),
        "study"=>data(s),"selected"=>data(selected),
        "evidence"=>[Dict("family"=>string(r.case.family),"level"=>r.case.level,
            "controls"=>data(controls(r.case)),"directory"=>abspath(dirname(r.history))) for r in results])
end
function publish_selection(path,s,selected,results)
    payload=selection_payload(s,selected,results)
    payload["selection_sha256"]=digest(payload)
    save_toml(path,payload)
end
function decode_aerodynamic_study(d)
    return merge(named(d),(;physical=named(d["physical"]),core_mode=Symbol(d["core_mode"]),
        nominal=named(d["nominal"]),sweeps=named(d["sweeps"]),families=Symbol.(d["families"]),
        absolute_tolerances=(CT=d["absolute_tolerances"]["CT"],CQ=d["absolute_tolerances"]["CQ"]),
        periodic_tolerances=(CT=d["periodic_tolerances"]["CT"],CQ=d["periodic_tolerances"]["CQ"])))
end
function load_selection(path)
    isfile(path) || error("Run chang_convergence_aerodynamic.jl first; selection missing: $path")
    payload=TOML.parsefile(path)
    get(payload,"schema",0)==1 && get(payload,"status","")=="verified" || error("Aerodynamic selection is not verified")
    seal=pop!(payload,"selection_sha256","")
    digest(payload)==seal || error("Aerodynamic selection file changed; regenerate it from the study")
    payload["source_sha256"]==SOURCE() || error("Numerical source changed since aerodynamic verification")
    s=decode_aerodynamic_study(payload["study"]); planned=validate_aerodynamic(s)
    results=NamedTuple[]
    for e in payload["evidence"]
        c=case(Symbol(e["family"]),e["level"],named(e["controls"]))
        options=aerodynamic_options(s.physical,s.core_mode,c,s)
        r=read_aerodynamic(e["directory"],options,s)
        push!(results,merge(r,(;case=c,reused=true,elapsed_s=0.)))
    end
    ordinary=filter(r->r.case.family!=:combined,results)
    sort([(string(r.case.family),r.case.level,key(r.case)) for r in ordinary])==
        sort([(string(c.family),c.level,key(c)) for c in planned]) || error("Selection evidence does not match the complete aerodynamic plan")
    selected=selected_controls(annotate(ordinary,planned,s),s)
    selected==controls(named(payload["selected"])) || error("Selected controls differ from the aerodynamic evidence")
    confirmation=confirmation_cases(selected,s)
    combined=sort(filter(r->r.case.family==:combined,results);by=r->r.case.level)
    length(combined)==2 && key.(getproperty.(combined,:case))==key.(confirmation) &&
        all(r->r.status=="completed",combined) && agrees(combined[1],combined[2],s) ||
        error("Combined aerodynamic refinement is not verified")
    return (;physical=s.physical,core_mode=s.core_mode,selected,source_sha256=SOURCE(),selection_sha256=seal)
end
function run_aerodynamic(s)
    cases=validate_aerodynamic(s); directory=abspath(s.output_directory)
    selection_path=joinpath(directory,"aerodynamic_selection.toml")
    if !s.dry_run
        if isfile(selection_path)
            mv(selection_path,selection_path*".previous_$(time_ns())";force=false)
        end
        save_toml(selection_path,(;schema=1,status="incomplete"))
    end
    results=execute_cases(cases,s,s.physical,s.core_mode;
        runner=(c,folder)->aerodynamic_case(c,s,folder))
    s.dry_run && return (;directory,cases,results,selection_file=nothing)
    rows=annotate(results,cases,s)
    selection_file=nothing; message="Selection not requested"
    if s.confirm_selection
        try
            selected=selected_controls(rows,s)
            confirmation=confirmation_cases(selected,s)
            for c in confirmation; aerodynamic_options(s.physical,s.core_mode,c,s); end
            # The combined cases are additional numerical evidence, never a
            # blind mixture of individually accepted controls.
            fingerprint=SOURCE()
            for c in confirmation
                println("Combined verification $(c.label): $(controls(c))")
                start=time(); r=aerodynamic_case(c,s,joinpath(directory,c.label))
                SOURCE()==fingerprint || error("Source changed during combined verification")
                push!(results,merge(r,(;case=c,elapsed_s=time()-start)))
            end
            append!(cases,confirmation); rows=annotate(results,cases,s)
            combined=filter(r->r.case.family==:combined,rows)
            all(r->r.accepted,combined) || error("Combined settings still change CT/CQ; refine the study before damping analysis")
            selection_file=publish_selection(selection_path,s,selected,results)
            # Re-read and independently recompute every acceptance decision.
            load_selection(selection_file)
            message="Verified selection: $selection_file"
        catch exception
            exception isa InterruptException && rethrow()
            message=sprint(showerror,exception)
            save_toml(selection_path,(;schema=1,status="not_converged",reason=message))
            selection_file=nothing
        end
    end
    write_tables(directory,rows,cases,s.physical,s.core_mode)
    write(joinpath(directory,"report.md"),"# Aerodynamic convergence\n\n$message\n\nCT/CQ means and phase waveforms control acceptance; CL remains in raw histories for inspection. The selected core is a numerical sensitivity choice, not a calibrated vortex radius. No retrimming is performed.\n")
    if s.make_plots
        @eval using Plots
        Base.invokelatest(plot_results,directory,rows)
    end
    println(message); println("Results: $directory")
    return (;directory,cases,results=rows,selection_file)
end
