module PazyLinearAeroelastic

let examples = normpath(joinpath(@__DIR__,".."))
    examples in LOAD_PATH || pushfirst!(LOAD_PATH,examples)
end

using LinearAlgebra
import AeroBeams
import WingPropellerUVLM as UVLM
include(joinpath(@__DIR__,"PazyWingUVLMCoupling.jl"))
include(joinpath(@__DIR__,"PazyLinearUVLM.jl"))

export linearize_pazy_aeroelastic,sweep_pazy_linear_flutter,
    save_pazy_coupled_modes,save_pazy_flutter_sweep,
    save_pazy_flutter_boundaries

function pazy_span_weights(pazy,spanwise_stations)
    node_count = pazy.n_elements + 1
    weights = zeros(spanwise_stations,node_count)
    for station in 1:spanwise_stations
        eta = (station - 1) / (spanwise_stations - 1)
        left = min(searchsortedlast(pazy.beamNodes,eta),pazy.n_elements)
        fraction = (eta - pazy.beamNodes[left]) /
                   (pazy.beamNodes[left+1] - pazy.beamNodes[left])
        weights[station,left] = 1 - fraction
        weights[station,left+1] = fraction
    end
    return weights
end

function pazy_grid(model,pazy,weights,chordwise_panels)
    return first(wing_geometry(model,pazy.chord,pazy.sparFraction,
        chordwise_panels,weights;sweepAngle=pazy.parameters.sweepAngle))
end

function grid_direction(model,x0,mode,pazy,weights,chordwise_panels,
    base_grid,base_offsets,step)

    displacement = Vector{Array{Float64,3}}(undef,2)
    vortex_offsets = Vector{Array{Float64,3}}(undef,2)
    for (part,component) in enumerate((real,imag))
        direction = component.(mode)
        if all(iszero,direction)
            displacement[part] = zeros(size(base_grid))
            vortex_offsets[part] = zeros(size(base_offsets))
            continue
        end
        geometries = map((-1.0,1.0)) do sign
            perturbed_model = deepcopy(model)
            AeroBeams.update_states!(perturbed_model,x0 .+ (sign * step) .* direction,
                1.0,1,1.0,0.0,Inf,false)
            wing_geometry(perturbed_model,pazy.chord,pazy.sparFraction,
                chordwise_panels,weights;sweepAngle=pazy.parameters.sweepAngle)
        end
        displacement[part] = (geometries[2][1] .- geometries[1][1]) ./ (2 * step)
        vortex_offsets[part] = (geometries[2][3] .- geometries[1][3]) ./ (2 * step)
    end
    return displacement,vortex_offsets
end

function residual_load_gain(problem,node_count)
    residual_reference = copy(AeroBeams.f!(problem))
    gain = zeros(length(residual_reference),6 * node_count)
    load_bcs = problem.model.BCs[2:end]
    length(load_bcs) == node_count || throw(DimensionMismatch(
        "Expected one six-component aerodynamic load BC per Pazy node"))

    try
        for node in 1:node_count,component in 1:6
            bc = load_bcs[node]
            previous_value = bc.state.currentValue[component]
            bc.state.currentValue[component] = previous_value + 1.0
            gain[:,6 * (node - 1) + component] .=
                AeroBeams.f!(problem) .- residual_reference
            bc.state.currentValue[component] = previous_value
        end
    finally
        AeroBeams.f!(problem)
    end
    return gain
end

function force_vector_to_nodal(force_vector,grid,vortex_offsets,weights)
    chordwise_vertices,spanwise_vertices = size(grid,2),size(grid,3)
    force_count = 3 * chordwise_vertices * spanwise_vertices
    node_count = size(weights,2)
    length(force_vector) == force_count || throw(DimensionMismatch(
        "UVLM force vector does not match the Pazy grid"))
    forces = Matrix{UVLM.SVector{3,Float64}}(undef,
        chordwise_vertices,spanwise_vertices)
    aerodynamic_loads = (;forces=[forces])
    span_map = PazySpanMap(weights)
    nodal_loads = zeros(6,node_count)

    for point in eachindex(forces)
        first_component = 3 * (point - 1) + 1
        forces[point] = UVLM.SVector{3,Float64}(force_vector[first_component],
            force_vector[first_component+1],force_vector[first_component+2])
    end
    beam_loads!(nodal_loads,aerodynamic_loads,vortex_offsets,span_map)
    return vec(nodal_loads)
end

function force_transfer_matrix(grid,vortex_offsets,weights)
    force_count = 3 * size(grid,2) * size(grid,3)
    transfer = zeros(6 * size(weights,2),force_count)
    unit_force = zeros(force_count)
    for column in 1:force_count
        unit_force[column] = 1.0
        transfer[:,column] .= force_vector_to_nodal(unit_force,grid,vortex_offsets,weights)
        unit_force[column] = 0.0
    end
    return transfer
end

function select_structural_modes(eigen_problem,number_of_modes,minimum_frequency_hz)
    eigenvalues = eigen_problem.eigensolution.λ
    positive = findall(λ -> isfinite(real(λ)) && isfinite(imag(λ)) &&
        imag(λ) > 2π * minimum_frequency_hz,eigenvalues)
    sort!(positive;by=index -> imag(eigenvalues[index]))
    length(positive) >= number_of_modes || throw(ArgumentError(
        "Requested $number_of_modes structural modes but only $(length(positive)) positive-frequency modes were found"))
    positive = positive[1:number_of_modes]
    negative = Int[]
    for index in positive
        partner = argmin(abs.(eigenvalues .- conj(eigenvalues[index])))
        partner != index && abs(eigenvalues[partner] - conj(eigenvalues[index])) <=
            1e-6 * max(abs(eigenvalues[index]),1.0) || throw(ArgumentError(
                "Structural eigenvalue $(eigenvalues[index]) has no conjugate partner"))
        push!(negative,partner)
    end
    return positive,negative
end

function linearize_pazy_aeroelastic(reference;
    number_of_modes=6,minimum_frequency_hz=1e-4,
    finite_difference_step=1e-6,kinematic_step=1e-5,
    steady_tolerance=1e-3,eigenmethod=:fulldecomp)

    reference.converged || throw(ArgumentError("Static aeroelastic reference did not converge"))
    number_of_modes isa Integer && number_of_modes > 0 || throw(ArgumentError(
        "number_of_modes must be a positive integer"))
    hasproperty(reference,:warm_start) || throw(ArgumentError(
        "Static reference is missing its Pazy geometry metadata"))

    aerodynamic = reference.aerodynamic_problem
    pazy = reference.warm_start.pazy
    base_system = aerodynamic.state.system
    chordwise_panels = size(base_system.grids[1],2) - 1
    spanwise_panels = size(base_system.grids[1],3) - 1
    weights = pazy_span_weights(pazy,spanwise_panels + 1)

    structural_problem = reference.structural_problem
    model = deepcopy(structural_problem.model)
    initial_state = copy(structural_problem.systemArrays.x)
    eigen_problem = AeroBeams.create_EigenProblem(
        model=model,eigenmethod=eigenmethod,
        nev=model.assemblyData.systemOrder,x0=initial_state)
    AeroBeams.solve!(eigen_problem)

    positive,negative = select_structural_modes(eigen_problem,number_of_modes,
        minimum_frequency_hz)
    eigenvalues = eigen_problem.eigensolution.λ
    descriptor_jacobian = Matrix(eigen_problem.systemArrays.jacobian)
    descriptor_inertia = Matrix(eigen_problem.systemArrays.inertia)
    left_spectrum = eigen(adjoint(descriptor_jacobian),-adjoint(descriptor_inertia))
    position_indices = reduce(vcat,
        vcat(element.assemblyIndices.DOF_u,element.assemblyIndices.DOF_p,
            element.assemblyIndices.DOF_F,element.assemblyIndices.DOF_M)
        for element in model.assemblyData.elements)

    positive_values = eigenvalues[positive]
    positive_vectors = Matrix{ComplexF64}(undef,length(initial_state),number_of_modes)
    positive_left = similar(positive_vectors)
    directions = Vector{Vector{Array{Float64,3}}}()
    base_grid,_,base_offsets = wing_geometry(model,pazy.chord,pazy.sparFraction,
        chordwise_panels,weights;sweepAngle=pazy.parameters.sweepAngle)
    maximum_absolute_difference(base_grid,base_system.grids[1]) / pazy.chord <= 1e-8 ||
        throw(ArgumentError("Structural and aerodynamic reference geometries differ"))
    offset_directions = Array{Float64,3}[]

    for (mode,index) in enumerate(positive)
        vector = eigen_problem.eigensolution.V[:,index]
        phase = -angle(sum(view(vector,position_indices).^2)) / 2
        vector .*= cis(phase)
        left_index = argmin(abs.(left_spectrum.values .- conj(eigenvalues[index])))
        left_vector = left_spectrum.vectors[:,left_index] .* cis(phase)
        overlap = dot(left_vector,descriptor_inertia * vector)
        abs(overlap) > eps() || throw(ArgumentError(
            "Could not biorthogonally normalize structural mode $(mode)"))
        left_vector ./= conj(overlap)
        positive_vectors[:,mode] .= vector
        positive_left[:,mode] .= left_vector

        mode_directions,mode_offset_directions = grid_direction(model,initial_state,
            vector,pazy,weights,chordwise_panels,base_grid,base_offsets,kinematic_step)
        push!(directions,[mode_directions[1]])
        push!(directions,[mode_directions[2]])
        push!(offset_directions,mode_offset_directions[1])
        push!(offset_directions,mode_offset_directions[2])
    end

    uv_structural = PazyLinearUVLM.linearize_uvlm_step(aerodynamic;
        displacement_directions=directions,velocity_directions=directions,
        finite_difference_step,steady_tolerance)
    mode_count = number_of_modes
    state_count = 2 * mode_count
    input_directions = 2 * mode_count
    motion_gain = zeros(ComplexF64,2 * input_directions,state_count)
    for mode in 1:mode_count
        positive_state = mode
        negative_state = mode_count + mode
        real_direction = 2 * mode - 1
        imaginary_direction = 2 * mode
        motion_gain[real_direction,positive_state] = 1
        motion_gain[imaginary_direction,positive_state] = im
        motion_gain[real_direction,negative_state] = 1
        motion_gain[imaginary_direction,negative_state] = -im
        velocity_offset = input_directions
        λ = positive_values[mode]
        motion_gain[velocity_offset + real_direction,positive_state] = λ
        motion_gain[velocity_offset + imaginary_direction,positive_state] = im * λ
        motion_gain[velocity_offset + real_direction,negative_state] = conj(λ)
        motion_gain[velocity_offset + imaginary_direction,negative_state] = -im * conj(λ)
    end

    base_grid,_,vortex_offsets = wing_geometry(model,pazy.chord,pazy.sparFraction,
        chordwise_panels,weights;sweepAngle=pazy.parameters.sweepAngle)
    force_transfer = force_transfer_matrix(base_grid,vortex_offsets,weights)
    load_gain = residual_load_gain(eigen_problem,pazy.n_elements + 1)
    modal_force_gain = Matrix{ComplexF64}(undef,state_count,size(uv_structural.C,1))
    geometric_force_gain = Matrix{ComplexF64}(undef,state_count,input_directions)
    time_step = uv_structural.time_step
    sampled_eigenvalues = vcat(exp.(positive_values .* time_step),
        exp.(conj.(positive_values) .* time_step))
    reference_load_transfer = Matrix{Float64}(undef,
        6 * (pazy.n_elements + 1),input_directions)
    for direction in 1:input_directions
        nodal_load_derivative = force_vector_to_nodal(uv_structural.force_reference,
            base_grid,offset_directions[direction],weights)
        reshape(nodal_load_derivative,6,:)[1:3,:] .= 0.0
        reference_load_transfer[:,direction] .= nodal_load_derivative
    end

    for mode in 1:mode_count
        continuous_gain = -vec(adjoint(positive_left[:,mode]) * load_gain * force_transfer)
        geometric_gain = -vec(adjoint(positive_left[:,mode]) * load_gain *
            reference_load_transfer)
        λ = positive_values[mode]
        integration_gain = abs(λ) > eps() ? (exp(λ * time_step) - 1) / λ : time_step
        modal_force_gain[mode,:] .= integration_gain .* continuous_gain
        modal_force_gain[mode_count + mode,:] .=
            conj(integration_gain .* continuous_gain)
        geometric_force_gain[mode,:] .= integration_gain .* geometric_gain
        geometric_force_gain[mode_count + mode,:] .= conj(integration_gain .* geometric_gain)
    end

    position_gain = motion_gain[1:input_directions,:]
    coupled_matrix = vcat(
        hcat(uv_structural.A,uv_structural.B * motion_gain),
        hcat(modal_force_gain * uv_structural.C,
            Diagonal(sampled_eigenvalues) + modal_force_gain * uv_structural.D * motion_gain +
            geometric_force_gain * position_gain)
    )
    coupled_eigensolution = eigen(coupled_matrix)
    discrete_eigenvalues = coupled_eigensolution.values
    continuous_eigenvalues = log.(discrete_eigenvalues) ./ time_step
    frequencies_hz = abs.(imag.(continuous_eigenvalues)) ./ (2π)
    growth_rates = real.(continuous_eigenvalues)
    damping_ratios = map(eachindex(continuous_eigenvalues)) do index
        ω = imag(continuous_eigenvalues[index])
        -100 * growth_rates[index] / hypot(growth_rates[index],ω)
    end

    return (;analysis=:coupled_frozen_wake,eigen_problem,uvlm=uv_structural,
        structural_eigenvalues=positive_values,discrete_eigenvalues,
        coupled_eigenvectors=coupled_eigensolution.vectors,
        continuous_eigenvalues,frequencies_hz,growth_rates_per_s=growth_rates,
        damping_ratios_percent=damping_ratios,geometric_force_gain,
        coupled_state_matrix=coupled_matrix)
end

function match_eigenmodes(previous_values,previous_vectors,current_values,current_vectors)
    length(previous_values) == length(current_values) || throw(DimensionMismatch(
        "Coupled state dimension changed during the speed sweep"))
    size(previous_vectors,1) == size(current_vectors,1) || throw(DimensionMismatch(
        "Coupled eigenvector dimensions changed during the speed sweep"))
    mode_count = length(previous_values)
    assigned = Vector{Int}(undef,mode_count)
    available = trues(mode_count)

    for previous_index in 1:mode_count
        previous_vector = view(previous_vectors,:,previous_index)
        previous_norm = max(norm(previous_vector),eps())
        best_index,best_cost = 0,Inf
        for current_index in 1:mode_count
            available[current_index] || continue
            current_vector = view(current_vectors,:,current_index)
            value_distance = abs(previous_values[previous_index] - current_values[current_index]) /
                max(abs(previous_values[previous_index]),abs(current_values[current_index]),1.0)
            overlap = abs(dot(previous_vector,current_vector)) /
                (previous_norm * max(norm(current_vector),eps()))
            cost = value_distance + 1 - min(overlap,1.0)
            if cost < best_cost
                best_index,best_cost = current_index,cost
            end
        end
        assigned[previous_index] = best_index
        available[best_index] = false
    end
    return assigned
end

function growth_rate_crossings(speeds,tracked_eigenvalues)
    size(tracked_eigenvalues,1) == length(speeds) || throw(DimensionMismatch(
        "Tracked eigenvalue rows must match the airspeed count"))
    boundaries = NamedTuple[]
    for mode in axes(tracked_eigenvalues,2),case_index in 1:length(speeds)-1
        lower_rate = real(tracked_eigenvalues[case_index,mode])
        upper_rate = real(tracked_eigenvalues[case_index+1,mode])
        lower_rate * upper_rate < 0 || continue
        lower_speed,upper_speed = speeds[case_index],speeds[case_index+1]
        crossing_speed = lower_speed - lower_rate *
            (upper_speed - lower_speed) / (upper_rate - lower_rate)
        push!(boundaries,(;mode,transition=lower_rate < 0 ? :onset : :restabilization,
            lower_speed_mps=lower_speed,upper_speed_mps=upper_speed,
            estimate_mps=crossing_speed))
    end
    return boundaries
end

"""Linearize converged references at increasing airspeeds and continue modes."""
function sweep_pazy_linear_flutter(airspeeds,reference_at_speed;
    number_of_modes=6,minimum_frequency_hz=1e-4,
    finite_difference_step=1e-6,kinematic_step=1e-5,
    steady_tolerance=1e-3,eigenmethod=:fulldecomp)

    speeds = Float64.(airspeeds)
    length(speeds) >= 2 && all(isfinite,speeds) && all(>(0),speeds) &&
        all(>(0),diff(speeds)) || throw(ArgumentError(
            "Airspeeds must be finite, positive, strictly increasing, and contain at least two values"))

    cases = map(speeds) do speed
        reference = reference_at_speed(speed)
        isapprox(reference.airspeed,speed;rtol=1e-10,atol=1e-10) ||
            throw(ArgumentError("Reference speed $(reference.airspeed) does not match requested speed $speed"))
        result = linearize_pazy_aeroelastic(reference;
            number_of_modes,minimum_frequency_hz,finite_difference_step,
            kinematic_step,steady_tolerance,eigenmethod)
        (;speed_mps=speed,result)
    end

    mode_count = length(first(cases).result.continuous_eigenvalues)
    tracked = Matrix{ComplexF64}(undef,length(cases),mode_count)
    tracked[1,:] .= first(cases).result.continuous_eigenvalues
    previous_vectors = first(cases).result.coupled_eigenvectors
    for case_index in 2:length(cases)
        current = cases[case_index].result
        order = match_eigenmodes(tracked[case_index-1,:],previous_vectors,
            current.continuous_eigenvalues,current.coupled_eigenvectors)
        tracked[case_index,:] .= current.continuous_eigenvalues[order]
        previous_vectors = current.coupled_eigenvectors[:,order]
    end

    boundaries = growth_rate_crossings(speeds,tracked)

    return (;airspeeds=speeds,cases,tracked_continuous_eigenvalues=tracked,boundaries)
end

function save_pazy_flutter_sweep(path::AbstractString,sweep)
    mkpath(dirname(path))
    open(path,"w") do io
        println(io,"speed_mps,mode,growth_rate_per_s,frequency_hz")
        for speed_index in eachindex(sweep.airspeeds),
            mode in axes(sweep.tracked_continuous_eigenvalues,2)
            λ = sweep.tracked_continuous_eigenvalues[speed_index,mode]
            println(io,join((string(sweep.airspeeds[speed_index]),string(mode),
                string(real(λ)),string(abs(imag(λ))/(2π))),','))
        end
    end
    return path
end

function save_pazy_flutter_boundaries(path::AbstractString,sweep)
    mkpath(dirname(path))
    open(path,"w") do io
        println(io,"mode,transition,lower_speed_mps,upper_speed_mps,estimate_mps")
        for boundary in sweep.boundaries
            println(io,join((string(boundary.mode),string(boundary.transition),
                string(boundary.lower_speed_mps),string(boundary.upper_speed_mps),
                string(boundary.estimate_mps)),','))
        end
    end
    return path
end

function save_pazy_coupled_modes(path::AbstractString,result)
    mkpath(dirname(path))
    open(path,"w") do io
        println(io,"mode,eigenvalue_real,eigenvalue_imag,frequency_hz,growth_rate_per_s,damping_ratio_percent")
        for index in eachindex(result.discrete_eigenvalues)
            λ = result.continuous_eigenvalues[index]
            println(io,join((string(index),string(real(λ)),string(imag(λ)),
                string(result.frequencies_hz[index]),string(result.growth_rates_per_s[index]),
                string(result.damping_ratios_percent[index])),','))
        end
    end
    return path
end

end