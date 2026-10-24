module PazyBifurcation

include(joinpath(@__DIR__, "..", "bifurcation", "TimeDomainBifurcation.jl"))
using .TimeDomainBifurcation
using Distributed, Dates
include(joinpath(@__DIR__, "PazyWingUVLMStatic.jl"))
export run_pazy_bifurcation, plot_bifurcation

const RESPONSE_DEFINITIONS = (
    tip_out_of_plane_m=(unit="m", label="Tip out-of-plane displacement", atol=1e-6),
    tip_out_of_plane_percent=(unit="% semispan", label="Tip out-of-plane displacement", atol=1e-4),
    tip_in_plane_m=(unit="m", label="Tip in-plane displacement (+downstream)", atol=1e-6),
    tip_twist_deg=(unit="deg", label="Tip twist", atol=1e-3),
    root_strain_ue=(unit="microstrain", label="Near-root longitudinal strain", atol=0.1),
    tip_geometric_incidence_deg=(unit="deg", label="Root incidence + elastic tip twist", atol=1e-3))

function model_responses(model, angle, span, gauge_y, gauge_z)
    tip = model.elements[end].nodalStates
    root = model.elements[1].compStates
    twist = wingtip_twist_degrees(model)
    return (; tip_out_of_plane_m=-tip.u_n2[1],
        tip_out_of_plane_percent=-100tip.u_n2[1]/span,
        tip_in_plane_m=-tip.u_n2[2], tip_twist_deg=twist,
        root_strain_ue=1e6*(root.γ[1] + gauge_z*root.κ[2] - gauge_y*root.κ[3]),
        tip_geometric_incidence_deg=rad2deg(angle)+twist)
end

function simulate_pazy_speed(speed, options)
    (; simulation_time, initialization, density, angle_of_attack, sideslip,
        chordwise_panels, spanwise_panels, time_step_chords, wake_length_chords,
        maximum_wake_rows, core_radius, wake_shedding_fraction, aerodynamic_interaction,
        symmetric_wing, pulse_amplitude, pulse_start, pulse_duration,
        gauge_y, gauge_z, dynamic_options, equilibrium_options) = options
    _, span, chord, _ = AeroBeams.geometrical_properties_Pazy()
    # Align the physical end time exactly, while respecting the requested U*dt/c.
    steps = ceil(Int, simulation_time/(time_step_chords*chord/speed))
    dt = simulation_time/steps
    steps >= 2 || throw(ArgumentError("Each simulation needs at least two timesteps"))
    if !iszero(pulse_amplitude)
        pulse_duration >= 8dt*(1-1e-12) || throw(ArgumentError("Resolve the tip pulse with at least eight timesteps"))
    end
    scaled_dt = dt*speed/chord
    wake_rows = isnothing(maximum_wake_rows) ?
        ceil(Int, wake_length_chords*chord/(speed*dt)) : maximum_wake_rows
    common = (; airspeed=speed, density, angle_of_attack,
        chordwise_panels, spanwise_panels, symmetric_wing,
        time_step_chords=scaled_dt, maximum_wake_rows=wake_rows,
        core_radius, wake_shedding_fraction, aerodynamic_interaction)
    # Both calls construct fresh problems. Never seed from another speed's result.
    equilibrium = initialization == :equilibrium ?
        run_pazy_wing_uvlm_static(; equilibrium_options..., common...) : nothing
    if !isnothing(equilibrium)
        equilibrium.converged || error("Static initialization failed at U=$speed: $(equilibrium.termination_reason)")
    end
    static = isnothing(equilibrium) ? nothing :
        model_responses(equilibrium.structural_problem.model, angle_of_attack, span, gauge_y, gauge_z)
    in_plane = Float64[isnothing(static) ? 0.0 : static.tip_in_plane_m]
    callback(i, time, structure, aerodynamic) = push!(in_plane,
        -structure.model.elements[end].nodalStates.u_n2[2])
    result = run_pazy_wing_uvlm(; dynamic_options..., common..., sideslip,
        initial_equilibrium=equilibrium, duration=simulation_time,
        initial_airspeed_fraction=1.0, airspeed_ramp_duration=0.0,
        settling_time=pulse_start, perturbation_amplitude=pulse_amplitude,
        perturbation_duration=pulse_duration, accepted_step_callback=callback,
        save_uvlm_history=false, uvlm_save_frequency=1, animation_frames=2,
        save_animation_history=false, save_structural_history=false,
        save_aerodynamic_load_history=false)
    all(==(speed), result.airspeed_history) || error("Bifurcation simulation must have constant airspeed")
    responses = (; tip_out_of_plane_m=result.tip_out_of_plane,
        tip_out_of_plane_percent=100 .* result.tip_out_of_plane ./ span,
        tip_in_plane_m=in_plane, tip_twist_deg=result.tip_twist_degrees,
        root_strain_ue=1e6 .* (result.root_axial_strain .+
            gauge_z .* result.root_curvature_2 .- gauge_y .* result.root_curvature_3),
        tip_geometric_incidence_deg=rad2deg(angle_of_attack) .+ result.tip_twist_degrees)
    return (; time=result.time, responses,
        static_equilibrium=isnothing(static) ? (;) : static)
end

"""Independent fixed-speed Pazy simulations using the existing nonlinear beam/UVLM solver.

Defaults: 35:1:52 m/s, 10 s each, statistics over t >= 8 s. Each speed starts
from its own converged static equilibrium and developed wake, with zero structural
velocity, followed by the same existing smooth tip-force pulse. `:undeformed`
instead uses the existing unloaded startup and an empty wake. No continuation.
The pulse coefficient is in N (peak approximately 0.6495*coefficient).

The Pazy model has no propeller DOFs. The reusable TimeDomainBifurcation module
accepts arbitrary scalar histories for other models, including pitch/yaw.
`worker_ids` optionally selects existing Distributed workers; the module is
loaded on those processes before pmap, and only compact summary rows return.
"""
function run_pazy_bifurcation(; enabled=true, U_min=35.0, U_max=52.0, dU=1.0,
    simulation_time=10.0, transient_time=8.0,
    responses=(:tip_out_of_plane_percent, :tip_in_plane_m, :tip_twist_deg, :root_strain_ue),
    initialization=:equilibrium, density=1.225, angle_of_attack=deg2rad(5.0), sideslip=0.0,
    chordwise_panels=10, spanwise_panels=20, time_step_chords=0.25,
    wake_length_chords=5.0, maximum_wake_rows=nothing,
    core_radius=1e-3, wake_shedding_fraction=0.1,
    aerodynamic_interaction=true, symmetric_wing=iszero(sideslip),
    pulse_amplitude=0.1, pulse_start=0.1, pulse_duration=0.04,
    gauge_y=0.0, gauge_z=2.25e-3/2,
    dynamic_options=(;), equilibrium_options=(;),
    response_definitions=RESPONSE_DEFINITIONS,
    check_convergence=true, rtol=0.02, convergence_cycles=3,
    save_time_histories=false, plot_results=true, plot_amplitude=true,
    zoom_limits=nothing, worker_ids=Int[],
    output_directory=joinpath(@__DIR__, "output", "bifurcation",
        Dates.format(now(), "yyyymmdd_HHMMSSsss")))
    enabled || return nothing
    bifurcation_speeds(U_min, U_max, dU)
    all(isfinite, (simulation_time, transient_time, density, angle_of_attack, sideslip,
        time_step_chords, wake_length_chords, pulse_amplitude, pulse_start,
        pulse_duration, gauge_y, gauge_z)) && 0 <= transient_time < simulation_time ||
        throw(ArgumentError("Require finite settings and 0 <= transient_time < simulation_time"))
    density > 0 && time_step_chords > 0 && wake_length_chords > 0 ||
        throw(ArgumentError("Density, timestep and wake length must be positive"))
    all(n -> n isa Integer && n > 0, (chordwise_panels, spanwise_panels)) ||
        throw(ArgumentError("Panel counts must be positive integers"))
    isnothing(maximum_wake_rows) || (maximum_wake_rows isa Integer && maximum_wake_rows > 0) ||
        throw(ArgumentError("maximum_wake_rows must be nothing or a positive integer"))
    initialization in (:equilibrium, :undeformed) || throw(ArgumentError("Unknown initialization"))
    initialization == :equilibrium && !iszero(sideslip) &&
        throw(ArgumentError("The existing static Pazy solver supports zero sideslip only"))
    !symmetric_wing || iszero(sideslip) || throw(ArgumentError("Symmetry requires zero sideslip"))
    0 <= pulse_start < simulation_time && pulse_duration > 0 ||
        throw(ArgumentError("Invalid pulse timing"))
    iszero(pulse_amplitude) || pulse_start+pulse_duration <= transient_time ||
        throw(ArgumentError("The excitation pulse must finish before the observation window"))
    all(n -> n in keys(RESPONSE_DEFINITIONS), responses) ||
        throw(ArgumentError("Select available Pazy responses: $(keys(RESPONSE_DEFINITIONS))"))
    UVLM.create_UVLMSolver(; coreRadius=core_radius,
        wakeSheddingFraction=wake_shedding_fraction, interaction=aerodynamic_interaction)
    # Only numerical controls may be overridden; state, speed and history options
    # remain owned by this independent-run driver.
    dynamic_defaults = (; newton_maximum_iterations=20, newton_absolute_tolerance=1e-5,
        newton_relative_tolerance=1e-5, newton_display_iterations=false,
        newton_always_update_jacobian=true, coupling_scheme=:strong,
        coupling_maximum_iterations=20, coupling_relaxation=0.5,
        coupling_geometry_tolerance=1e-4, coupling_load_tolerance=1e-4,
        coupling_display_iterations=false, progress_frequency=100)
    static_defaults = (; maximum_aerodynamic_time=0.20, minimum_aerodynamic_time=0.05,
        aerodynamic_window_time=0.01, aerodynamic_force_tolerance=1e-4,
        aerodynamic_circulation_tolerance=1e-4, maximum_static_iterations=100,
        static_relaxation=0.3, geometry_tolerance=1e-5, load_tolerance=1e-4,
        consecutive_equilibrium_iterations=2, newton_maximum_iterations=50,
        newton_absolute_tolerance=1e-8, newton_relative_tolerance=1e-8,
        newton_display_iterations=false, progress_frequency=100)
    for (given, defaults) in ((dynamic_options, dynamic_defaults), (equilibrium_options, static_defaults))
        all(k -> k in keys(defaults), keys(given)) ||
            throw(ArgumentError("Only numerical solver controls may be overridden"))
    end
    options = (; simulation_time, initialization, density, angle_of_attack, sideslip,
        chordwise_panels, spanwise_panels, time_step_chords, wake_length_chords,
        maximum_wake_rows, core_radius, wake_shedding_fraction, aerodynamic_interaction,
        symmetric_wing, pulse_amplitude, pulse_start, pulse_duration, gauge_y, gauge_z,
        dynamic_options=merge(dynamic_defaults, dynamic_options),
        equilibrium_options=merge(static_defaults, equilibrium_options))
    # Module loading precedes serialization of the model-specific closure.
    all(id -> id != myid() && id in workers(), worker_ids) ||
        throw(ArgumentError("Select active remote worker IDs"))
    if !isempty(worker_ids)
        Distributed.remotecall_eval(Main, worker_ids, quote
            if !isdefined(Main, :PazyBifurcation)
                Base.include(Main, $(joinpath(@__DIR__, "PazyBifurcation.jl")))
            end
        end)
    end
    simulate = speed -> simulate_pazy_speed(speed, options)
    return run_bifurcation(simulate; U_min, U_max, dU,
        definitions=response_definitions, responses, simulation_time, transient_time,
        check_convergence, rtol, convergence_cycles, save_time_histories,
        plot_results, plot_amplitude, zoom_limits, worker_ids, metadata=options, output_directory)
end

end # module
