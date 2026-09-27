# Diagnostic only: run an instrumented copy of the current coupling function.
# Neither solver nor the production runner is modified. Outputs go to a new folder.
# Usage: julia --project=lib/WingPropellerUVLM audit_pazy_failure.jl loose 3.0
# Modes: loose, strong, loose_freshjac (same case, every Newton Jacobian refreshed).
using LinearAlgebra, Serialization, Dates

const AUDIT_MODE = isempty(ARGS) ? "loose" : ARGS[1]
const AUDIT_DURATION = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 3.0
const AUDIT_DIR = joinpath(@__DIR__, "output", "failure_audit_" *
    Dates.format(now(), "yyyymmdd_HHMMSS") * "_" * AUDIT_MODE)
mkpath(AUDIT_DIR)
const AUDIT_IO = open(joinpath(AUDIT_DIR, "history.tsv"), "w")
const AUDIT_PROBES = open(joinpath(AUDIT_DIR, "residual_checks.txt"), "w")
const AUDIT_BEFORE = Ref{Any}(nothing)
const AUDIT_CONTEXT = Ref{Any}(nothing)
const AUDIT_SOLVES = Ref(0)
const AUDIT_START = time()

# Exact copies with diagnostic hooks; saved structural/animation histories are
# omitted to keep the failure checkpoints small. Those saves do not drive motion.
source_path = joinpath(@__DIR__, "PazyWingUVLMCoupling.jl")
source = read(source_path, String)
source = replace(source, "function run_pazy_wing_uvlm(" => "function audit_run(")
source = replace(source, "        AeroBeams.get_equivalent_states_rates!(structure)" =>
    "        audit_begin!(structure, aerodynamic, grid, previous_grid, vortex_offsets, previous_offsets, weights, i)\n        AeroBeams.get_equivalent_states_rates!(structure)")
source = replace(source, "AeroBeams.solve_time_step!(structure)" => "audit_solve!(structure)")
source = replace(source, "        AeroBeams.save_time_step_data!(structure, time[i])" =>
    "        audit_record!(structure, aerodynamic, candidate_grid, accepted_aerodynamic_loads, geometry_residual, coupling_iterations)")
source = replace(source, "save_wing_frame!(aerodynamic, surface_history, wake_history, animation_time)" => "nothing")
source = replace(source, "            UVLM.rollback_time_step!(aerodynamic)" =>
    "            audit_failure!(structure, aerodynamic)\n            UVLM.rollback_time_step!(aerodynamic)")
Base.include_string(Main, source, source_path)

function audit_begin!(structure, aerodynamic, grid, previous_grid, offsets, previous_offsets, weights, i)
    AUDIT_CONTEXT[] = (; structure, aerodynamic, grid, previous_grid, offsets, previous_offsets, weights, i)
    AUDIT_SOLVES[] = 0
    # Keep the accepted state before a late failure, including the wake at n.
    if structure.timeNow >= 2.7 || i % 500 == 0
        AUDIT_BEFORE[] = deepcopy(AUDIT_CONTEXT[])
    end
end

function residual_blocks(p)
    elements = p.model.elements
    ids(names) = unique(vcat([getproperty(e, name) for e in elements for name in names]...))
    return (; force=norm(p.residual[ids((:eqs_Fu1, :eqs_Fu2))]),
        moment=norm(p.residual[ids((:eqs_Fp1, :eqs_Fp2))]),
        compatibility=norm(p.residual[ids((:eqs_FF1, :eqs_FF2, :eqs_FM1, :eqs_FM2))]),
        velocity=norm(p.residual[ids((:eqs_FV, :eqs_FΩ))]))
end

function audit_solve!(p)
    AUDIT_SOLVES[] += 1
    sample = AUDIT_SOLVES[] == 1 && (p.indexEndTimeStep % 200 == 0)
    if sample
        probe = deepcopy(p)
        AeroBeams.assemble_system_arrays!(probe)
        println(AUDIT_PROBES, "t=", p.timeNow, " before=", norm(probe.residual), " ", residual_blocks(probe))
    end
    AeroBeams.solve_time_step!(p)
    if sample || !p.systemSolver.convergedFinalSolution
        probe = deepcopy(p)
        AeroBeams.assemble_system_arrays!(probe)
        println(AUDIT_PROBES, "t=", p.timeNow, " converged=", p.systemSolver.convergedFinalSolution,
            " reported=", norm(p.residual), " reassembled=", norm(probe.residual), " ", residual_blocks(probe))
        flush(AUDIT_PROBES)
    end
end

panel_area(p) = norm(cross(p.rbr-p.rtl, p.rtr-p.rbl))/2

function audit_record!(structure, aerodynamic, structural_grid, loads, mismatch, iterations)
    sys = aerodynamic.state.system
    wake = @view sys.wakes[1][1:sys.nwake[1], :]
    # Compare the displacement increment actually used in UVLM with the same
    # increment of the structural surface, both expressed on the vortex lattice.
    ctx = AUDIT_CONTEXT[]
    _, _, previous_panels = UVLM.grid_to_surface_panels(ctx.grid)
    _, _, actual_panels = UVLM.grid_to_surface_panels(structural_grid)
    aero_points = UVLM.imperial_nodal_positions(sys.surfaces[1])
    old_aero_points = UVLM.imperial_nodal_positions(sys.previous_surfaces[1])
    struct_points = UVLM.imperial_nodal_positions(actual_panels)
    old_struct_points = UVLM.imperial_nodal_positions(previous_panels)
    fa = UVLM.dimensional_loads(sys).forces[1]
    work_aero = sum(dot(fa[k], aero_points[k]-old_aero_points[k]) for k in eachindex(fa))
    work_struct_surface = sum(dot(fa[k], struct_points[k]-old_struct_points[k]) for k in eachindex(fa))
    # This is a displacement-based interface-work diagnostic, not exact beam
    # energy or a proof of fluid/structure energy conservation.
    vals = (structure.timeNow, -structure.model.elements[end].nodalStates.u_n2[1],
        wingtip_twist_degrees(structure.model), mismatch, iterations,
        maximum(norm, eachcol(loads[1:3,:])), maximum(norm, eachcol(loads[4:6,:])),
        maximum(abs, sys.Γ), maximum(abs, sys.dΓdt), maximum(norm, sys.Vcp[1]),
        maximum(norm, sys.V[1]), minimum(panel_area, sys.surfaces[1]), minimum(panel_area, wake),
        maximum(norm(e.states.V) for e in structure.model.elements),
        maximum(norm(e.states.Ω) for e in structure.model.elements),
        norm(structure.residual), work_aero, work_struct_surface, work_struct_surface-work_aero)
    println(AUDIT_IO, join(vals, '\t'))
    if structure.indexEndTimeStep % 200 == 0
        flush(AUDIT_IO)
        println("AUDIT t=", round(structure.timeNow; digits=4), " mismatch/c=", round(mismatch; sigdigits=4),
            " wall=", round(time()-AUDIT_START; digits=1), " s")
    end
end

function audit_failure!(structure, aerodynamic)
    serialize(joinpath(AUDIT_DIR, "failure.jls"),
        (; before=AUDIT_BEFORE[], failed=structure, aerodynamic, context=AUDIT_CONTEXT[]))
    println("Failure checkpoint saved in ", AUDIT_DIR)
end

function main()
    @assert AUDIT_MODE in ("loose", "strong", "loose_freshjac")
    println(AUDIT_IO, "time\ttip_bending_m\ttip_twist_deg\tgeometry_mismatch_over_c\tfsi_iterations\tmax_nodal_force_N\tmax_nodal_moment_Nm\tmax_gamma\tmax_gammadot\tmax_surface_motion_speed\tmax_wake_speed\tmin_surface_area\tmin_wake_area\tmax_beam_speed\tmax_beam_angular_speed\tnr_reported_residual\taero_increment_work_J\tstruct_surface_increment_work_J\twork_difference_J")
    settings = (; airspeed=50.0, density=1.225, angle_of_attack=deg2rad(3.0), sideslip=0.0,
        initial_airspeed_fraction=1.0, airspeed_ramp_duration=0.0,
        duration=AUDIT_DURATION, settling_time=min(1.0, AUDIT_DURATION/2),
        chordwise_panels=10, spanwise_panels=20, time_step_chords=0.25,
        symmetric_wing=true, maximum_wake_rows=40, core_radius=1e-3,
        wake_shedding_fraction=0.1, aerodynamic_interaction=true,
        save_uvlm_history=false, uvlm_save_frequency=1,
        newton_maximum_iterations=20, newton_absolute_tolerance=1e-6,
        newton_relative_tolerance=1e-6, newton_display_iterations=false,
        newton_always_update_jacobian=AUDIT_MODE=="loose_freshjac",
        perturbation_amplitude=0.0, perturbation_duration=0.1,
        animation_frames=2, animation_time_step=10.0, progress_frequency=500,
        coupling_scheme=AUDIT_MODE=="strong" ? :strong : :loose,
        coupling_maximum_iterations=20, coupling_relaxation=0.5,
        coupling_geometry_tolerance=1e-5, coupling_load_tolerance=1e-5,
        coupling_display_iterations=false)
    open(joinpath(AUDIT_DIR, "settings.txt"), "w") do io
        println(io, settings)
        println(io, "Julia ", VERSION, "; threads=", Threads.nthreads(), "; BLAS threads=", BLAS.get_num_threads())
    end
    try
        result = audit_run(; settings...)
        serialize(joinpath(AUDIT_DIR, "result.jls"), result)
        println("AUDIT completed: ", AUDIT_DIR)
    catch err
        open(joinpath(AUDIT_DIR, "error.txt"), "w") do io
            showerror(io, err, catch_backtrace())
        end
        showerror(stdout, err)
        println()
    finally
        close(AUDIT_IO)
        close(AUDIT_PROBES)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
