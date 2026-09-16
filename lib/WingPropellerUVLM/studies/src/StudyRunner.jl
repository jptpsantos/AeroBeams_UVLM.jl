function operating_selection(selection,s)
    speed=get(s,:speed_mps,nothing)
    isnothing(speed) && return selection
    speed isa Real && !(speed isa Bool) && isfinite(speed) && speed>0 ||
        throw(ArgumentError("Aeroelastic speed_mps must be finite and positive, or nothing"))
    p=selection.physical
    rpm=p.rpm*(speed/p.speed_mps)
    isfinite(rpm) && rpm>0 || throw(ArgumentError("Scaled aeroelastic RPM must be finite and positive"))
    # The verification seal belongs to the original operating point only.
    return (;physical=merge(p,(;speed_mps=Float64(speed),rpm)),
        core_mode=selection.core_mode,selected=selection.selected,
        aerodynamic_reference=selection)
end

"""Run damping convergence with an optional speed override after verifying the source selection."""
function run_aeroelastic(s)
    check_loaded_source()
    reference=Convergence.load_selection(s.selection_file)
    selection=operating_selection(reference,s)
    println("Aeroelastic operating point: $(selection.physical.speed_mps) m/s, $(selection.physical.rpm) RPM")
    Convergence.load_model()
    result=Base.invokelatest(Convergence.run_aeroelastic_loaded,MovingBlockStudy(s),selection)
    description="Aerodynamic reference: $(reference.physical.speed_mps) m/s, $(reference.physical.rpm) RPM.\n\nAeroelastic operating point: $(selection.physical.speed_mps) m/s, $(selection.physical.rpm) RPM.\n\n"
    if selection.physical.speed_mps!=reference.physical.speed_mps
        description*="The aerodynamic selection was verified at the reference speed; this study checks aeroelastic sensitivity at the overridden speed with proportional RPM.\n\n"
    end
    open(joinpath(result.directory,"report.md"),"a") do io
        print(io,"\n",description)
    end
    return merge(result,(;aerodynamic_reference=reference))
end
