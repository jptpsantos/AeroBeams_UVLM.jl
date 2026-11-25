if !isdefined(@__MODULE__,:PazyFlutterSweep)
    include(joinpath(@__DIR__,"PazyFlutterSweep.jl"))
end
using .PazyFlutterSweep: run_pazy_flutter_sweep
using LinearAlgebra: BLAS

function pazy_flutter_configuration(; pilot_only=false,angle_deg=nothing)
    commissioningSlice = (angle_deg=5.,pilot_speed=45.,
        speeds=(40.,43.,44.,45.,49.,50.,50.5,52.),lco_seed=45.,lco_limits=(35.,53.))
    # Within each incidence, these are search speeds in m/s. The historical
    # results seed 0–7 degrees. The 8–10 degree grids are exploratory.
    slices = (
        (angle_deg=0.,pilot_speed=70.,speeds=(55.,60.,63.,64.,65.,67.,70.,75.,80.,90.,100.,110.,120.),
            lco_seed=70.,lco_limits=(58.,90.)),
        (angle_deg=1.,pilot_speed=67.,speeds=(53.,57.,59.,60.,63.,67.,72.,75.,76.,80.,84.),
            lco_seed=67.,lco_limits=(53.,83.)),
        (angle_deg=2.,pilot_speed=58.,speeds=(47.,51.,53.,54.,56.,58.,62.,65.,66.,70.,73.),
            lco_seed=58.,lco_limits=(45.,73.)),
        (angle_deg=3.,pilot_speed=53.,speeds=(45.,48.,50.,53.,56.,58.,60.),
            lco_seed=53.,lco_limits=(39.,62.)),
        (angle_deg=4.,pilot_speed=49.,speeds=(41.,44.,45.,46.,47.,49.,51.,53.,54.,57.),
            lco_seed=49.,lco_limits=(38.,58.)),
        commissioningSlice,
        (angle_deg=6.,pilot_speed=44.,speeds=(36.,39.,40.,41.,42.,44.,45.,46.,47.,49.,51.),
            lco_seed=44.,lco_limits=(33.,51.)),
        (angle_deg=7.,pilot_speed=42.,speeds=(35.,38.,39.,40.,42.,44.,45.,48.),
            lco_seed=42.,lco_limits=(34.,48.)),
        (angle_deg=8.,pilot_speed=40.,speeds=(30.,34.,36.,38.,39.,40.,41.,42.,44.,46.,49.),
            lco_seed=40.,lco_limits=(29.,49.)),
        (angle_deg=9.,pilot_speed=39.,speeds=(28.,32.,34.,36.,38.,39.,40.,42.,44.,48.),
            lco_seed=39.,lco_limits=(28.,48.)),
        (angle_deg=10.,pilot_speed=38.,speeds=(26.,30.,32.,34.,36.,37.,38.,40.,42.,44.,48.),
            lco_seed=38.,lco_limits=(26.,48.)),
    )
    if !isnothing(angle_deg)
        slices = filter(slice -> slice.angle_deg == angle_deg,slices)
        isempty(slices) && throw(ArgumentError("No configured slice at $angle_deg deg"))
    elseif pilot_only
        slices = (commissioningSlice,)
    end
    return (;
        slices,
        model_options=(; mesh=:refined,rootTorsionCorr=false),
        # Use the 0.2 ms setting that trimmed at 5 deg and 45 m/s for the
        # initial scan. The edge check halves this step. Keep 0.1 ms in the pilot.
        chordwise_panels=10,spanwise_panels=20,coarse_dt=0.0004,
        fine_dt=pilot_only ? 0.0001 : 0.0002,
        wake_length_chords=6.5,wake_factor=2.0,
        duration=2.0,pulse_forces=(0.05,0.10),
        growth_options=(; windows=((0.50,1.20),(1.20,1.90)),frequency_band=(22.0,43.0)),
        speed_tolerance=0.5,lco_step=2.0,plateau_time=10.0,
        include_entry_paths=true,refine_lco=true,save_histories=true
    )
end

function main_pazy_flutter_sweep(args=ARGS)
    valid = ("--pilot","--equilibrium","--lco","--all","--check-config","--help")
    remaining = collect(args)
    angleIndex = findfirst(==("--angle"),remaining)
    angle = nothing
    if !isnothing(angleIndex)
        angleIndex < length(remaining) || throw(ArgumentError("--angle requires an incidence in degrees"))
        angle = tryparse(Float64,remaining[angleIndex+1])
        !isnothing(angle) && isfinite(angle) || throw(ArgumentError("Invalid --angle value"))
        deleteat!(remaining,angleIndex:angleIndex+1)
    end
    outputIndex = findfirst(==("--output-directory"),remaining)
    outputDirectory = nothing
    if !isnothing(outputIndex)
        outputIndex < length(remaining) || throw(ArgumentError("--output-directory requires a path"))
        path = remaining[outputIndex+1]
        !isempty(path) && !startswith(path,"--") || throw(ArgumentError("--output-directory requires a path"))
        outputDirectory = abspath(path)
        deleteat!(remaining,outputIndex:outputIndex+1)
    end
    length(remaining) <= 1 && all(in(valid),remaining) ||
        throw(ArgumentError("Use one of $(join(valid,", ")), optionally with --angle DEG and --output-directory PATH"))
    option = isempty(remaining) ? "--equilibrium" : only(remaining)
    if option == "--help"
        println("Pazy flutter sweep. Edit pazy_flutter_configuration() to select the model and speed grids.")
        println("--equilibrium (default) runs fixed-speed pulse cases and refines flutter brackets.")
        println("--lco runs stepped LCO paths. --all runs both sweeps.")
        println("--pilot runs only the timestep and wake diagnostic at 5 deg by default.")
        println("--angle DEG selects one configured incidence. --check-config displays the schedule.")
        println("--output-directory PATH saves the run to a new directory at PATH.")
        println("Use julia --threads=10 to start Julia with 10 threads.")
        return nothing
    end
    configuration = pazy_flutter_configuration(; pilot_only=option == "--pilot",angle_deg=angle)
    if option == "--check-config"
        check = run_pazy_flutter_sweep(; configuration...,stage=:equilibrium,dry_run=true)
        println("$(length(check.slices)) incidences, $(sum(s -> length(s.speeds),check.slices)) fixed-speed cases before refinement.")
        println("Sweep dt=$(check.options.fine_dt) s, retained wake=$(check.options.wake_length_chords) chords.")
        println("Julia threads: ",Threads.nthreads())
        println("Output directory: ",something(outputDirectory,joinpath(@__DIR__,"output","flutter_sweep","run_*")))
        for (slice,paths) in zip(check.slices,check.paths)
            println("Root incidence $(slice.angle_deg) deg, pilot at $(slice.pilot_speed) m/s.")
            println("Equilibrium speeds [m/s]: ",slice.speeds)
            for path in paths
                sequence = path.direction == :downsweep ? reverse(path.speeds) : path.speeds
                println(path.name," [m/s]: ",sequence)
            end
        end
        return check
    end

    # Small dense factorizations run inside the threaded UVLM calculation.
    previousBLASThreads = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    try
        flutter_results = run_pazy_flutter_sweep(; configuration...,stage=Symbol(option[3:end]),
            output_directory=outputDirectory)
        if option == "--pilot"
            passed = count(p -> p.status == :passed,flutter_results.pilots)
            println("Numerical pilot passed at $passed of $(length(flutter_results.pilots)) incidences.")
        end
        cases = flutter_results.equilibrium_cases
        if !isempty(cases)
            accepted = count(c -> c.status in (:growth,:decay),cases)
            failed = count(c -> c.status == :failed,cases)
            println("Fixed-speed evaluations: $accepted accepted, $(length(cases)-accepted-failed) unresolved, $failed failed.")
        end
        println("Boundary intervals recorded: ",length(flutter_results.boundaries))
        println("Outputs: ",flutter_results.output_directory)
        return flutter_results
    finally
        BLAS.set_num_threads(previousBLASThreads)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    flutter_results = main_pazy_flutter_sweep()
end
