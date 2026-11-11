# Bounded synthetic-wake benchmark. Run with --threads=4 and the package project.
# Writes a report into a fresh directory; never launches a full simulation.
module ChangRuntimeBenchmark

using LinearAlgebra,StaticArrays,Printf,Statistics
import WingPropellerUVLM as W
const EXAMPLE = normpath(joinpath(@__DIR__,".."))
include(joinpath(EXAMPLE,"src","ChangAeroelastic.jl"))
include(joinpath(EXAMPLE,"src","ChangRuntimeOptimizations.jl"))
const C = ChangAeroelastic
const F = ChangRuntimeOptimizations

function populate_wakes!(model,workspace)
    C.update_aero_geometry_for_state!(model,workspace,zeros(model.structural.ndof_free),0.0)
    p = model.parameters
    for isurf in eachindex(workspace.system.surfaces)
        panels = workspace.system.surfaces[isurf]
        for j in axes(workspace.system.wakes[isurf],2),i in axes(workspace.system.wakes[isurf],1)
            left,right = W.bottom_left(panels[end,j]),W.bottom_right(panels[end,j])
            function aged(point,age)
                translation = SVector(p.Vinf * age,0.0,0.0)
                isurf == 1 && return point + translation
                ip = div(isurf - 2,p.Nb_prop) + 1
                hub = workspace.T_hub_A_current[ip]
                return hub + W.RotationMatrix(p.Ω * age,1) * (point - hub) + translation
            end
            a,b = (i - 1) * p.dt[1],i * p.dt[1]
            workspace.system.wakes[isurf][i,j] = W.WakePanel(
                aged(left,a),aged(right,a),aged(left,b),aged(right,b),
                W.get_core_size(panels[end,j]),1.0 + 0.001i + 0.01j)
        end
    end
    for isurf in eachindex(workspace.system.surfaces)
        surface = workspace.system.surfaces[isurf]
        workspace.system.wake_shedding_locations[isurf] .=
            vcat([W.bottom_left(surface[end,j]) for j in axes(surface,2)],
                [W.bottom_right(surface[end,end])])
    end
    workspace.system.Γ .= 1.0
end

function measure_commit(f,input,workspace,dt,active_rows)
    system = deepcopy(input)
    stats = @timed f(system,system.freestream[],dt;
        repeated_points=workspace.repeated_points,nwake=active_rows,
        interaction_id=workspace.surface_interaction_id,
        interaction=false)
    return (; seconds=stats.time,bytes=stats.bytes,system)
end

function main()
    BLAS.set_num_threads(1)
    config = C.load_chang_configuration(; env=Dict{String,String}())
    model = C.build_chang_model(config)
    workspace = F.build_fast_workspace(model)
    populate_wakes!(model,workspace)
    system,dt = workspace.system,model.parameters.dt[1]
    println("Warming up original and threaded wake commits; threads=$(Threads.nthreads())")
    flush(stdout)
    for f in (W.advance_wake!,F.advance_fast_wake!)
        measure_commit(f,system,workspace,dt,ones(Int,workspace.nsurf))
    end
    serial,threaded = Float64[],Float64[]
    for repetition = 1:3
        a = measure_commit(W.advance_wake!,system,workspace,dt,workspace.nwake)
        b = measure_commit(F.advance_fast_wake!,system,workspace,dt,workspace.nwake)
        @assert a.system.V == b.system.V
        @assert a.system.wakes == b.system.wakes
        push!(serial,a.seconds)
        push!(threaded,b.seconds)
        @printf("Repetition %d: original %.4f s, threaded %.4f s, speedup %.2fx; exact equality verified\n",
            repetition,a.seconds,b.seconds,a.seconds / b.seconds)
        flush(stdout)
    end
    base = mkpath(joinpath(EXAMPLE,"output","validation"))
    directory = mktempdir(base;
        prefix="fast_runtime_benchmark_",cleanup=false)
    report = joinpath(directory,"report.md")
    open(report,"w") do io
        println(io,"# Optional Chang wake benchmark\n")
        println(io,"Julia $VERSION; $(Threads.nthreads()) Julia threads; BLAS threads: 1.\n")
        println(io,"Synthetic fully populated wakes; current mesh and wake capacities; interaction disabled. ",
            "Compilation warmup and copies excluded. Other simulations may affect wall time. ",
            "This is a wake-commit benchmark, not an end-to-end simulation or remaining-time estimate.\n")
        println(io,"Bound panels: $(length(system.Γ)); wake panels: $(sum(length,system.wakes)); ",
            "wake rows: $(workspace.nwake).\n")
        println(io,"| Repetition | Original serial (s) | Threaded (s) |\n|---|---:|---:|")
        for i in eachindex(serial)
            @printf(io,"| %d | %.6f | %.6f |\n",i,serial[i],threaded[i])
        end
        @printf(io,
            "\nMedian wake speedup: **%.2fx**. Velocities and shed/convected wake panels matched exactly in all repetitions.\n",
            median(serial) / median(threaded))
        bytesPerStep = length(system.Γ) * sizeof(W.SurfacePanel{Float64})
        @printf(io,
            "\nDisabling full surface history avoids %.2f GiB of retained panel data over the current %d-step case. ",
            bytesPerStep * length(model.parameters.dt) / 2.0^30,length(model.parameters.dt))
        println(io,"Animation frames and the structural response history remain available independently.")
    end
    println("Benchmark report: $report")
end

main()
end # module
