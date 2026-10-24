module WingPropellerBifurcation

# Reuse the general case solver with its existing structural integrator choices.
const UVLM_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
UVLM_ROOT in LOAD_PATH || pushfirst!(LOAD_PATH, UVLM_ROOT)
include(joinpath(@__DIR__, "src", "WingPropellerAeroelastic.jl"))
include(joinpath(@__DIR__, "..", "bifurcation", "TimeDomainBifurcation.jl"))
using .WingPropellerAeroelastic, .TimeDomainBifurcation
using WingPropellerUVLM: full_beam_aerodynamic_kinematics
using Distributed, Dates
export run_wing_propeller_bifurcation, plot_bifurcation

function response_definitions(case)
    pairs = Pair{Symbol,Any}[
        :tip_out_of_plane_m => (unit="m", label="Tip out-of-plane displacement (+up)", atol=1e-6),
        :tip_in_plane_m => (unit="m", label="Tip in-plane displacement (+chord)", atol=1e-6),
        :tip_twist_rad => (unit="rad", label="Tip twist about span", atol=1e-5)]
    for i in eachindex(case.propellers), direction in (:pitch, :yaw)
        push!(pairs, Symbol("propeller_$(i)_$(direction)_rad") =>
            (unit="rad", label="Propeller $i $direction", atol=1e-5))
    end
    return (; pairs...)
end

function simulate_wing_propeller_speed(speed, base_case, options)
    case = deepcopy(base_case)
    operating_condition = merge(case.operating_condition, (; freestream_speed=speed))
    case = merge(case, (; operating_condition))
    # Reassembly recalculates constant-J angular speeds at this U, together with
    # their gyroscopic terms. Fixed-RPM propellers retain the case's speed law.
    structural = build_structural_model(case)
    nominal_dt = isnothing(options.time_step_s) ?
        deg2rad(options.azimuth_step_deg)/maximum(abs, structural.angular_speeds) : options.time_step_s
    isfinite(nominal_dt) && nominal_dt > 0 || throw(ArgumentError("Invalid physical timestep"))
    steps = max(2, ceil(Int, options.simulation_time/nominal_dt))
    time = collect(0:steps) .* (options.simulation_time/steps)
    result = run_uvlm_time_domain_analysis(case; time, time_steps=diff(time),
        options.solver_options...)
    names = keys(response_definitions(case))
    responses = NamedTuple{names}(Tuple(zeros(length(time)) for _ in names))
    for (j,q) in enumerate(result.displacement_history)
        k = full_beam_aerodynamic_kinematics(result.structural, q)
        responses.tip_out_of_plane_m[j] = last(k.u_z)
        responses.tip_in_plane_m[j] = last(k.u_x)
        responses.tip_twist_rad[j] = last(k.theta_y)
        for i in eachindex(case.propellers), direction in (:pitch, :yaw)
            field = Symbol("propeller_", direction)
            name = Symbol("propeller_$(i)_$(direction)_rad")
            getproperty(responses, name)[j] = getproperty(k, field)[i]
        end
    end
    return (; time, responses, angular_speeds=structural.angular_speeds)
end

"""Independent runs of an existing Xu/Chang/Bohnisch case with pitch/yaw outputs.

The base case retains its physical properties, speed law and structural solver
options (Newmark or generalized-alpha, existing coupling). Supply the unchanged
UVLM runner keywords in `solver_options`, including copied initial displacement
or velocity and an external load callback when desired. Each speed constructs a
new model and starts with a new aerodynamic system and empty wake. The existing
solver retains displacement/velocity/acceleration histories during one run;
these are released after extracting scalar responses, and are not saved by default.
"""
function run_wing_propeller_bifurcation(case; enabled=true,
    U_min=60., U_max=110., dU=1., simulation_time=10., transient_time=8.,
    responses=keys(response_definitions(case)), time_step_s=nothing,
    azimuth_step_deg=5., solver_options=(;), worker_ids=Int[],
    output_directory=joinpath(@__DIR__, "output", "bifurcation",
        Dates.format(now(), "yyyymmdd_HHMMSSsss")), analysis_options...)
    enabled || return nothing
    isfinite(azimuth_step_deg) && azimuth_step_deg > 0 || throw(ArgumentError("azimuth_step_deg must be positive"))
    isnothing(time_step_s) || (isfinite(time_step_s) && time_step_s > 0) ||
        throw(ArgumentError("time_step_s must be nothing or positive"))
    any(k -> k in (:time, :time_steps, :aerodynamic_load, :begin_step, :commit_step), keys(solver_options)) &&
        throw(ArgumentError("Time and wake transactions are owned by the existing UVLM runner"))
    options = (; simulation_time, time_step_s, azimuth_step_deg, solver_options)
    all(id -> id != myid() && id in workers(), worker_ids) ||
        throw(ArgumentError("Select active remote worker IDs"))
    if !isempty(worker_ids)
        Distributed.remotecall_eval(Main, worker_ids, quote
            if !isdefined(Main, :WingPropellerBifurcation)
                Base.include(Main, $(joinpath(@__DIR__, "WingPropellerBifurcation.jl")))
            end
        end)
    end
    # Each task gets independent options as well as a fresh model and wake.
    simulate = speed -> simulate_wing_propeller_speed(speed, case, deepcopy(options))
    return run_bifurcation(simulate; U_min, U_max, dU, simulation_time, transient_time,
        responses, definitions=response_definitions(case), worker_ids, output_directory,
        metadata=(; case, options), analysis_options...)
end

end # module
