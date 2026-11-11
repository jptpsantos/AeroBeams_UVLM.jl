# Isolated-propeller trim solver shared by the windmilling and powered studies.

using LinearAlgebra
using Printf
using StaticArrays
using Statistics

using WingPropellerUVLM:
    Freestream,
    Reference,
    RotationMatrix,
    System,
    copy_surfaces_to_previous!,
    generate_propeller_blades_grid,
    get_nodal_properties_chang,
    get_twist_deg_interp,
    grid_to_surface_panels,
    imperial_nodal_forces,
    imperial_nodal_positions,
    legacy_imperial_segment_forces!,
    near_field_forces!,
    propagate_system!,
    repeated_trailing_edge_points

"""Numerical and physical controls for isolated Chang propeller trims."""
Base.@kwdef struct ChangWindmillingTrimOptions
    flow_speed_mps::Float64 = 65.0
    air_density_kgpm3::Float64 = 1.225
    angle_of_attack_deg::Float64 = 3.0
    sideslip_deg::Float64 = 0.0
    propeller_radius_m::Float64 = 1.15
    propeller_chord_m::Float64 = 0.197
    blade_count::Int = 4
    radial_panels::Int = 20
    chordwise_panels::Int = 20
    # Added uniformly to the complete Chang twist distribution. Leave at zero
    # to use the blade geometry defined in BladeGeometry.jl unchanged.
    collective_pitch_offset_deg::Float64 = 0.0
    azimuth_step_deg::Float64 = 5.0
    simulated_revolutions::Int = 4
    averaged_revolutions::Int = 1
    retained_wake_revolutions::Int = 3
    wake_relaxation::Float64 = 0.1
    core_radius_m::Union{Nothing,Float64} = nothing
    vortex_core_span_fraction::Float64 = 0.5
    vortex_core_chord_fraction::Float64 = 0.0
    near_field_force_model::Symbol = :imperial
    initial_rpm::Float64 = 1207.96
    rpm_bracket::NTuple{2,Float64} = (1200,1400)
    torque_tolerance_nm::Float64 = 0.25
    thrust_tolerance_n::Float64 = 1.0
    rpm_tolerance::Float64 = 0.25
    maximum_root_iterations::Int = 8
end

function validate_chang_windmilling_options(options::ChangWindmillingTrimOptions)
    options.flow_speed_mps > 0 || error("flow_speed_mps must be positive")
    options.air_density_kgpm3 > 0 || error("air_density_kgpm3 must be positive")
    options.propeller_radius_m > 0 || error("propeller_radius_m must be positive")
    options.propeller_chord_m > 0 || error("propeller_chord_m must be positive")
    options.blade_count > 0 || error("blade_count must be positive")
    options.radial_panels > 0 || error("radial_panels must be positive")
    options.chordwise_panels > 0 || error("chordwise_panels must be positive")
    isfinite(options.collective_pitch_offset_deg) ||
        error("collective_pitch_offset_deg must be finite")
    options.azimuth_step_deg > 0 || error("azimuth_step_deg must be positive")
    options.simulated_revolutions >= 2 || error("simulate at least two revolutions")
    1 <= options.averaged_revolutions < options.simulated_revolutions ||
        error("averaged_revolutions must be positive and smaller than simulated_revolutions")
    options.retained_wake_revolutions > 0 || error("retained_wake_revolutions must be positive")
    options.wake_relaxation >= 0 || error("wake_relaxation must be nonnegative")
    if isnothing(options.core_radius_m)
        factors = (options.vortex_core_span_fraction,options.vortex_core_chord_fraction)
        all(x -> isfinite(x) && x >= 0,factors) ||
            error("Vortex-core factors must be finite and nonnegative")
        max(factors...) > 0 || error("A positive finite core is required for free-wake self induction")
    else
        isfinite(options.core_radius_m) && options.core_radius_m > 0 ||
            error("core_radius_m must be finite and positive, or nothing to use core factors")
    end
    options.near_field_force_model in (:imperial,:legacy_imperial_segments) ||
        error(
            "near_field_force_model must be :imperial or " *
            ":legacy_imperial_segments"
        )
    options.initial_rpm > 0 || error("initial_rpm must be positive")
    0 < options.rpm_bracket[1] < options.rpm_bracket[2] ||
        error("rpm_bracket must contain two increasing positive values")
    options.torque_tolerance_nm > 0 || error("torque_tolerance_nm must be positive")
    options.thrust_tolerance_n > 0 || error("thrust_tolerance_n must be positive")
    options.rpm_tolerance > 0 || error("rpm_tolerance must be positive")
    options.maximum_root_iterations > 0 || error("maximum_root_iterations must be positive")

    steps_per_revolution = round(Int,360 / options.azimuth_step_deg)
    isapprox(
        steps_per_revolution * options.azimuth_step_deg,
        360.0;
        atol=100 * eps(Float64),
        rtol=0.0
    ) || error("azimuth_step_deg must divide 360 degrees exactly")
    return steps_per_revolution
end

"""Return the Chang blade angle at 75% radius, including any collective offset."""
chang_blade_angle_75_deg(options::ChangWindmillingTrimOptions) =
    get_twist_deg_interp(0.75) + options.collective_pitch_offset_deg

"""Return the near-field implementation selected by the trim options."""
function chang_trim_near_field_function(model::Symbol)
    model === :imperial && return near_field_forces!
    model === :legacy_imperial_segments && return legacy_imperial_segment_forces!
    error("Unsupported near-field force model: $model")
end

"""Core callback shared by initial and rotated trim geometry."""
function chang_trim_finite_core(options::ChangWindmillingTrimOptions)
    if isnothing(options.core_radius_m)
        return (chord,segmentWidth) -> max(
            options.vortex_core_span_fraction * segmentWidth,
            options.vortex_core_chord_fraction * chord
        )
    end
    radius = Float64(options.core_radius_m)
    return (chord,segmentWidth) -> radius
end

function chang_trim_core_description(options::ChangWindmillingTrimOptions)
    return isnothing(options.core_radius_m) ?
           "max($(options.vortex_core_span_fraction) Δs, $(options.vortex_core_chord_fraction) c)" :
           "$(options.core_radius_m) m (fixed)"
end

function rotate_chang_propeller_grids!(currentGrids,referenceGrids,omega,time)
    spin = RotationMatrix(-omega * time,1)
    for bladeIndex in eachindex(referenceGrids)
        referenceGrid = referenceGrids[bladeIndex]
        currentGrid = currentGrids[bladeIndex]
        for vertexIndex in CartesianIndices((size(referenceGrid,2),size(referenceGrid,3)))
            i,j = Tuple(vertexIndex)
            currentGrid[:,i,j] = spin * SVector{3}(referenceGrid[:,i,j])
        end
    end
    return currentGrids
end

"""Return total rotor thrust and aerodynamic shaft torque from nodal loads."""
function chang_propeller_shaft_loads(system; hub_position=SVector(0.0,0.0,0.0))
    nodal_forces = imperial_nodal_forces(system)
    nodalPositions = imperial_nodal_positions(system)
    totalForce = zero(hub_position)
    totalMoment = zero(hub_position)

    for surface_index in eachindex(nodal_forces)
        for vertexIndex in eachindex(nodal_forces[surface_index])
            force = nodal_forces[surface_index][vertexIndex]
            position = nodalPositions[surface_index][vertexIndex]
            totalForce += force
            totalMoment += cross(position - hub_position,force)
        end
    end

    # The UVLM body frame points +x downstream. Positive thrust is defined upstream.
    return (thrust_n=-totalForce[1],torque_nm=totalMoment[1])
end

"""
    simulate_chang_propeller_rpm(rpm, options; verbose=true)

March the isolated rigid Chang propeller with a free wake and average the
selected near-field shaft loads over complete final revolutions.
"""
function simulate_chang_propeller_rpm(
    rpm::Real,
    options::ChangWindmillingTrimOptions;
    verbose::Bool=true
)
    steps_per_revolution = validate_chang_windmilling_options(options)
    rpm > 0 || error("rpm must be positive")

    omega = Float64(rpm) * 2pi / 60
    dt = deg2rad(options.azimuth_step_deg) / omega
    total_steps = options.simulated_revolutions * steps_per_revolution
    maximum_wake_rows = options.retained_wake_revolutions * steps_per_revolution
    near_field_force_function = chang_trim_near_field_function(
        options.near_field_force_model
    )

    _,_,twists_at_nodes = get_nodal_properties_chang(options.radial_panels)
    blade_twists = deg2rad.(
        twists_at_nodes .- 90.0 .+ options.collective_pitch_offset_deg,
    )
    referenceGrids = generate_propeller_blades_grid(
        options.propeller_radius_m,
        options.propeller_chord_m,
        options.radial_panels,
        options.chordwise_panels,
        blade_twists,
        options.blade_count
    )
    currentGrids = deepcopy(referenceGrids)
    fcore = chang_trim_finite_core(options)
    initialSurfaces = [
        grid_to_surface_panels(grid; fcore=fcore)[3]
        for grid in currentGrids
    ]
    system = System(initialSurfaces; nw=fill(maximum_wake_rows,options.blade_count))
    system.reference[] = Reference(
        pi * options.propeller_radius_m^2,
        options.propeller_chord_m,
        2 * options.propeller_radius_m,
        SVector(0.0,0.0,0.0),
        options.flow_speed_mps,
        options.air_density_kgpm3
    )
    system.symmetric .= false
    system.surface_id .= 1:options.blade_count
    system.wake_finite_core .= true
    system.trailing_vortices .= false
    system.surfaces[:] = initialSurfaces
    system.previous_surfaces[:] = deepcopy(initialSurfaces)
    system.freestream[] = Freestream(
        options.flow_speed_mps,
        deg2rad(options.angle_of_attack_deg),
        deg2rad(options.sideslip_deg),
        SVector(0.0,0.0,0.0)
    )

    repeated_points = repeated_trailing_edge_points(system.surfaces)
    active_wake_rows = zeros(Int,options.blade_count)
    interaction_id = collect(1:options.blade_count)
    thrustHistory = Vector{Float64}(undef,total_steps)
    torqueHistory = Vector{Float64}(undef,total_steps)

    verbose && @printf(
        "  RPM=%8.3f: %d steps, %d retained wake rows, model=%s\n",
        rpm,
        total_steps,
        maximum_wake_rows,
        string(options.near_field_force_model),
    )

    for step = 1:total_steps
        copy_surfaces_to_previous!(system,options.blade_count)
        time = step * dt
        rotate_chang_propeller_grids!(currentGrids,referenceGrids,omega,time)
        for bladeIndex = 1:options.blade_count
            system.surfaces[bladeIndex] = grid_to_surface_panels(
                currentGrids[bladeIndex];
                fcore=fcore
            )[3]
        end

        propagate_system!(
            system,
            system.freestream[],
            dt;
            additional_velocity=nothing,
            repeated_points=repeated_points,
            nwake=active_wake_rows,
            eta=options.wake_relaxation,
            calculate_influence_matrix=true,
            near_field_analysis=true,
            near_field_force_function=near_field_force_function,
            derivatives=false,
            interaction_id=interaction_id,
            interaction=true,
            advance_wake=true
        )
        active_wake_rows .= min.(active_wake_rows .+ 1,maximum_wake_rows)

        loads = chang_propeller_shaft_loads(system)
        thrustHistory[step] = loads.thrust_n
        torqueHistory[step] = loads.torque_nm
        if verbose && step % steps_per_revolution == 0
            completedRevolution = step ÷ steps_per_revolution
            revolutionRange = (step-steps_per_revolution+1):step
            @printf(
                "      revolution %d/%d: mean Q=%+.4f N m\n",
                completedRevolution,
                options.simulated_revolutions,
                mean(torqueHistory[revolutionRange]),
            )
        end
    end

    revolutionMeanTorque = [
        mean(torqueHistory[((rev-1)*steps_per_revolution+1):(rev*steps_per_revolution)])
        for rev = 1:options.simulated_revolutions
    ]
    revolutionMeanThrust = [
        mean(thrustHistory[((rev-1)*steps_per_revolution+1):(rev*steps_per_revolution)])
        for rev = 1:options.simulated_revolutions
    ]
    averageSteps = options.averaged_revolutions * steps_per_revolution
    averageRange = (total_steps-averageSteps+1):total_steps
    meanTorque = mean(torqueHistory[averageRange])
    meanThrust = mean(thrustHistory[averageRange])
    torqueStandardDeviation = std(torqueHistory[averageRange])
    thrustStandardDeviation = std(thrustHistory[averageRange])
    nRevolutionsPerSecond = omega / (2pi)
    diameter = 2 * options.propeller_radius_m
    advance_ratio = options.flow_speed_mps / (nRevolutionsPerSecond * diameter)
    torque_coefficient = meanTorque / (
        options.air_density_kgpm3 * nRevolutionsPerSecond^2 * diameter^5
    )
    thrust_coefficient = meanThrust / (
        options.air_density_kgpm3 * nRevolutionsPerSecond^2 * diameter^4
    )
    periodicTorqueChange = revolutionMeanTorque[end] - revolutionMeanTorque[end-1]
    periodicThrustChange = revolutionMeanThrust[end] - revolutionMeanThrust[end-1]

    return (
        near_field_force_model=options.near_field_force_model,
        rpm=Float64(rpm),
        omega_radps=omega,
        advance_ratio=advance_ratio,
        inflow_ratio=advance_ratio / pi,
        mean_thrust_n=meanThrust,
        mean_torque_nm=meanTorque,
        thrust_coefficient=thrust_coefficient,
        torque_coefficient=torque_coefficient,
        torque_standard_deviation_nm=torqueStandardDeviation,
        thrust_standard_deviation_n=thrustStandardDeviation,
        periodic_torque_change_nm=periodicTorqueChange,
        periodic_thrust_change_n=periodicThrustChange,
        revolution_mean_torque_nm=revolutionMeanTorque,
        revolution_mean_thrust_n=revolutionMeanThrust,
        thrust_history_n=thrustHistory,
        torque_history_nm=torqueHistory,
        steps_per_revolution=steps_per_revolution,
        time_step_s=dt
    )
end

# Backward-compatible name used by earlier windmilling scripts.
simulate_chang_windmilling_rpm(args...; kwargs...) =
    simulate_chang_propeller_rpm(args...; kwargs...)

function _opposite_sign_or_zero(a,b)
    return iszero(a) || iszero(b) || signbit(a) != signbit(b)
end

"""Solve mean aerodynamic shaft torque = 0 with a safeguarded secant iteration."""
function trim_chang_windmilling_rpm(
    options::ChangWindmillingTrimOptions;
    verbose::Bool=true
)
    validate_chang_windmilling_options(options)
    evaluations = Dict{Float64,NamedTuple}()

    function evaluate(rpm)
        key = round(Float64(rpm); digits=8)
        result = get!(evaluations,key) do
            simulate_chang_propeller_rpm(key,options; verbose=verbose)
        end
        verbose && @printf(
            "    mean Q=%+10.4f N m, CQ=%+.6e, J=%.6f, rev drift=%+.4f N m\n",
            result.mean_torque_nm,
            result.torque_coefficient,
            result.advance_ratio,
            result.periodic_torque_change_nm,
        )
        return result
    end

    initialResult = evaluate(options.initial_rpm)
    if abs(initialResult.mean_torque_nm) <= options.torque_tolerance_nm
        return (trim=initialResult,initial=initialResult,evaluations=collect(values(evaluations)))
    end

    lowerRPM,upperRPM = options.rpm_bracket
    lowerResult = evaluate(lowerRPM)
    upperResult = evaluate(upperRPM)
    expansionCount = 0
    while !_opposite_sign_or_zero(lowerResult.mean_torque_nm,upperResult.mean_torque_nm)
        expansionCount += 1
        expansionCount <= 5 || error(
            "Unable to bracket zero torque after five RPM-bracket expansions"
        )
        lowerRPM = max(50.0,0.75 * lowerRPM)
        upperRPM = 1.25 * upperRPM
        lowerResult = evaluate(lowerRPM)
        upperResult = evaluate(upperRPM)
    end

    if abs(lowerResult.mean_torque_nm) <= options.torque_tolerance_nm
        orderedEvaluations = sort!(collect(values(evaluations)); by=result -> result.rpm)
        return (trim=lowerResult,initial=initialResult,evaluations=orderedEvaluations)
    elseif abs(upperResult.mean_torque_nm) <= options.torque_tolerance_nm
        orderedEvaluations = sort!(collect(values(evaluations)); by=result -> result.rpm)
        return (trim=upperResult,initial=initialResult,evaluations=orderedEvaluations)
    end

    bestResult = abs(lowerResult.mean_torque_nm) <= abs(upperResult.mean_torque_nm) ?
                 lowerResult : upperResult
    for _ = 1:options.maximum_root_iterations
        torqueSpan = upperResult.mean_torque_nm - lowerResult.mean_torque_nm
        trialRPM = upperRPM - upperResult.mean_torque_nm *
                              (upperRPM - lowerRPM) / torqueSpan
        bracketWidth = upperRPM - lowerRPM
        edgeMargin = 0.1 * bracketWidth
        if !isfinite(trialRPM) ||
           trialRPM <= lowerRPM + edgeMargin || trialRPM >= upperRPM - edgeMargin
            trialRPM = (lowerRPM + upperRPM) / 2
        end

        trialResult = evaluate(trialRPM)
        if abs(trialResult.mean_torque_nm) < abs(bestResult.mean_torque_nm)
            bestResult = trialResult
        end
        if abs(trialResult.mean_torque_nm) <= options.torque_tolerance_nm ||
           bracketWidth <= options.rpm_tolerance
            bestResult = trialResult
            break
        end

        if _opposite_sign_or_zero(lowerResult.mean_torque_nm,trialResult.mean_torque_nm)
            upperRPM,upperResult = trialRPM,trialResult
        else
            lowerRPM,lowerResult = trialRPM,trialResult
        end
    end

    abs(bestResult.mean_torque_nm) <= options.torque_tolerance_nm || @warn(
        "Windmilling trim stopped above the requested torque tolerance",
        residualTorqueNm = bestResult.mean_torque_nm,
        torque_tolerance_nm = options.torque_tolerance_nm,
    )
    orderedEvaluations = sort!(collect(values(evaluations)); by=result -> result.rpm)
    return (trim=bestResult,initial=initialResult,evaluations=orderedEvaluations)
end

"""
    trim_chang_thrusting_rpm(options; target_thrust_n, verbose=true)

Solve `mean_thrust_n = target_thrust_n` by varying RPM. Unlike windmilling
trim, the resulting aerodynamic shaft torque is generally nonzero and must be
balanced by the motor.
"""
function trim_chang_thrusting_rpm(
    options::ChangWindmillingTrimOptions;
    target_thrust_n::Real,
    verbose::Bool=true
)
    validate_chang_windmilling_options(options)
    isfinite(target_thrust_n) && target_thrust_n > 0 ||
        error("target_thrust_n must be finite and positive")
    targetThrust = Float64(target_thrust_n)
    evaluations = Dict{Float64,NamedTuple}()

    function evaluate(rpm)
        key = round(Float64(rpm); digits=8)
        result = get!(evaluations,key) do
            simulate_chang_propeller_rpm(key,options; verbose=verbose)
        end
        residual = result.mean_thrust_n - targetThrust
        verbose && @printf(
            "    mean T=%+10.4f N, residual=%+10.4f N, CT=%+.6e, J=%.6f, rev drift=%+.4f N\n",
            result.mean_thrust_n,
            residual,
            result.thrust_coefficient,
            result.advance_ratio,
            result.periodic_thrust_change_n,
        )
        return result
    end

    thrust_residual(result) = result.mean_thrust_n - targetThrust
    initialResult = evaluate(options.initial_rpm)
    if abs(thrust_residual(initialResult)) <= options.thrust_tolerance_n
        return (
            trim=initialResult,
            initial=initialResult,
            target_thrust_n=targetThrust,
            residual_thrust_n=thrust_residual(initialResult),
            evaluations=collect(values(evaluations))
        )
    end

    lowerRPM,upperRPM = options.rpm_bracket
    lowerResult = evaluate(lowerRPM)
    upperResult = evaluate(upperRPM)
    lowerResidual = thrust_residual(lowerResult)
    upperResidual = thrust_residual(upperResult)
    expansionCount = 0
    while !_opposite_sign_or_zero(lowerResidual,upperResidual)
        expansionCount += 1
        expansionCount <= 5 || error(
            "Unable to bracket the target thrust after five RPM-bracket expansions"
        )
        lowerRPM = max(50.0,0.75 * lowerRPM)
        upperRPM = 1.25 * upperRPM
        lowerResult = evaluate(lowerRPM)
        upperResult = evaluate(upperRPM)
        lowerResidual = thrust_residual(lowerResult)
        upperResidual = thrust_residual(upperResult)
    end

    if abs(lowerResidual) <= options.thrust_tolerance_n
        orderedEvaluations = sort!(collect(values(evaluations)); by=result -> result.rpm)
        return (
            trim=lowerResult,
            initial=initialResult,
            target_thrust_n=targetThrust,
            residual_thrust_n=lowerResidual,
            evaluations=orderedEvaluations
        )
    elseif abs(upperResidual) <= options.thrust_tolerance_n
        orderedEvaluations = sort!(collect(values(evaluations)); by=result -> result.rpm)
        return (
            trim=upperResult,
            initial=initialResult,
            target_thrust_n=targetThrust,
            residual_thrust_n=upperResidual,
            evaluations=orderedEvaluations
        )
    end

    bestResult = abs(lowerResidual) <= abs(upperResidual) ?
                 lowerResult : upperResult
    for _ = 1:options.maximum_root_iterations
        residualSpan = upperResidual - lowerResidual
        trialRPM = upperRPM - upperResidual *
                              (upperRPM - lowerRPM) / residualSpan
        bracketWidth = upperRPM - lowerRPM
        edgeMargin = 0.1 * bracketWidth
        if !isfinite(trialRPM) ||
           trialRPM <= lowerRPM + edgeMargin ||
           trialRPM >= upperRPM - edgeMargin
            trialRPM = (lowerRPM + upperRPM) / 2
        end

        trialResult = evaluate(trialRPM)
        trialResidual = thrust_residual(trialResult)
        if abs(trialResidual) < abs(thrust_residual(bestResult))
            bestResult = trialResult
        end
        if abs(trialResidual) <= options.thrust_tolerance_n ||
           bracketWidth <= options.rpm_tolerance
            bestResult = trialResult
            break
        end

        if _opposite_sign_or_zero(lowerResidual,trialResidual)
            upperRPM,upperResult,upperResidual =
                trialRPM,trialResult,trialResidual
        else
            lowerRPM,lowerResult,lowerResidual =
                trialRPM,trialResult,trialResidual
        end
    end

    finalResidual = thrust_residual(bestResult)
    abs(finalResidual) <= options.thrust_tolerance_n || @warn(
        "Thrusting trim stopped above the requested thrust tolerance",
        residual_thrust_n = finalResidual,
        thrust_tolerance_n = options.thrust_tolerance_n,
    )
    orderedEvaluations = sort!(collect(values(evaluations)); by=result -> result.rpm)
    return (
        trim=bestResult,
        initial=initialResult,
        target_thrust_n=targetThrust,
        residual_thrust_n=finalResidual,
        evaluations=orderedEvaluations
    )
end
