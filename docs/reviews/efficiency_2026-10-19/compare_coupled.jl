# Capture the original solver before changes, then compare the optimized solver.
# julia --startup-file=no --compiled-modules=existing --project=. <this file> capture
# julia --startup-file=no --compiled-modules=existing --project=. <this file> compare
using AeroBeams, LinearAlgebra, TOML, Test
include(joinpath(@__DIR__, "..", "..", "..", "lib", "WingPropellerUVLM",
    "examples", "pazy_aerobeams_coupled", "PazyWingUVLMCoupling.jl"))

function short_case(scheme; kwargs...)
    run_pazy_wing_uvlm(; airspeed=30.0, density=1.225, angle_of_attack=0.01,
        sideslip=0.0, initial_airspeed_fraction=1.0, airspeed_ramp_duration=0.0,
        duration=0.005, settling_time=0.0008, chordwise_panels=3, spanwise_panels=6,
        symmetric_wing=true, maximum_wake_rows=6, core_radius=1e-3,
        wake_shedding_fraction=0.1, aerodynamic_interaction=true,
        save_uvlm_history=false, uvlm_save_frequency=1,
        newton_maximum_iterations=20, newton_absolute_tolerance=1e-8,
        newton_relative_tolerance=1e-8, newton_display_iterations=false,
        newton_always_update_jacobian=true, perturbation_amplitude=0.05,
        perturbation_duration=0.002, animation_frames=2, progress_frequency=100,
        time_step_chords=0.25, save_animation_history=false, coupling_scheme=scheme,
        coupling_maximum_iterations=35, coupling_relaxation=0.5,
        coupling_geometry_tolerance=1e-7, coupling_load_tolerance=1e-5, kwargs...)
end

function record_result(result)
    sys = result.aerodynamic_problem.state.system
    rows = sys.nwake[1]
    wake = [component for panel in sys.wakes[1][1:rows,:]
        for point in (UVLM.top_left(panel),UVLM.top_right(panel),UVLM.bottom_left(panel),UVLM.bottom_right(panel))
        for component in point]
    Dict("tip_displacement"=>result.tip_out_of_plane,
        "tip_twist"=>result.tip_twist_degrees, "circulation"=>copy(sys.Γ),
        "state"=>copy(result.structural_problem.x), "wake"=>wake,
        "nodal_loads"=>vec(result.aerodynamic_nodal_load_history),
        "iterations"=>result.coupling_iterations,
        "time"=>result.time, "wake_rows"=>rows)
end

function main()
    path = joinpath(@__DIR__,"coupled_baseline.toml")
    mode = isempty(ARGS) ? "compare" : only(ARGS)
    mode in ("capture","compare") || error("Use capture or compare")
    results = Dict(string(scheme)=>short_case(scheme) for scheme in (:strong,:loose))
    records = Dict(scheme=>record_result(result) for (scheme,result) in results)
    if mode == "capture"
        isfile(path) && error("Baseline already exists; preserve the original reference")
        open(io->TOML.print(io,records),path,"w")
    else
        baseline = TOML.parsefile(path)
        @testset "Coupled solver equivalence to original source" begin
            for scheme in ("strong","loose")
                for field in ("tip_displacement","tip_twist","circulation","state","wake","nodal_loads")
                    @test isapprox(records[scheme][field],baseline[scheme][field];rtol=2e-8,atol=2e-10)
                end
                @test records[scheme]["iterations"] == baseline[scheme]["iterations"]
                @test records[scheme]["wake_rows"] == baseline[scheme]["wake_rows"]
                @test records[scheme]["time"] == baseline[scheme]["time"]
            end
        end
        @testset "Reduced history leaves the coupled response unchanged" begin
            for scheme in (:strong,:loose)
                reduced = short_case(scheme;save_structural_history=false,
                    save_aerodynamic_load_history=false)
                actual = record_result(reduced)
                expected = records[string(scheme)]
                for field in ("tip_displacement","tip_twist","circulation","state","wake")
                    @test isapprox(actual[field],expected[field];rtol=2e-8,atol=2e-10)
                end
                @test actual["iterations"] == expected["iterations"]
                @test isempty(reduced.structural_problem.xOverTime)
                @test isempty(reduced.structural_problem.nodalStatesOverTime)
                @test isempty(reduced.structural_problem.savedTimeVector)
                @test isempty(reduced.aerodynamic_nodal_load_history)
            end
        end
        @testset "Pazy load transfer preserves virtual work" begin
            weights = [1.0 0 0;0.3 0.7 0;0 0.4 0.6;0 0 1.0]
            offsets = reshape(sin.(1.0:27),3,3,3) .* 0.1
            forces = [UVLM.SVector(sin(k+j),cos(k-j),sin(2k-j)) for k in 1:3,j in 1:4]
            displacement = reshape(cos.(1.0:9),3,3)
            rotation = reshape(sin.(1.0:9),3,3)
            loads = beam_loads((;forces=[forces]),offsets,weights)
            structural_work = sum(loads[1:3,:] .* displacement)+sum(loads[4:6,:] .* rotation)
            aerodynamic_work = 0.0
            for j in 1:4,k in 1:3
                delta = zeros(3)
                for node in 1:3
                    delta .+= weights[j,node] .* (displacement[:,node]+cross(rotation[:,node],offsets[:,k,node]))
                end
                aerodynamic_work += dot(forces[k,j],A_TO_UVLM*delta)
            end
            @test structural_work ≈ aerodynamic_work atol=1e-13 rtol=1e-13
        end
    end
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
