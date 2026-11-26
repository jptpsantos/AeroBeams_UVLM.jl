module PazyLinearStructuralModes

let examples = normpath(joinpath(@__DIR__,".."))
    examples in LOAD_PATH || pushfirst!(LOAD_PATH,examples)
end

import AeroBeams
using DelimitedFiles

export solve_pazy_structural_modes,save_pazy_structural_modes

"""Solve AeroBeams modes at a converged static state with aerodynamic loads frozen."""
function solve_pazy_structural_modes(reference; number_of_modes::Integer=20,
    eigenmethod::Symbol=:krylovkit)

    number_of_modes > 0 || throw(ArgumentError("number_of_modes must be positive"))
    hasproperty(reference,:converged) ||
        throw(ArgumentError("reference must contain a convergence flag"))
    reference.converged || throw(ArgumentError("static aeroelastic reference did not converge"))
    hasproperty(reference,:structural_problem) ||
        throw(ArgumentError("reference must contain a structural_problem"))

    structural = reference.structural_problem
    model = deepcopy(structural.model)
    initial_state = copy(structural.systemArrays.x)
    isempty(initial_state) && throw(ArgumentError("static structural state is empty"))
    available_modes = model.assemblyData.systemOrder
    eigen_problem = AeroBeams.create_EigenProblem(
        model=model,
        eigenmethod=eigenmethod,
        nev=min(number_of_modes,available_modes),
        x0=initial_state
    )
    AeroBeams.solve!(eigen_problem)

    modes = map(eachindex(eigen_problem.eigensolution.λ)) do index
        λ = eigen_problem.eigensolution.λ[index]
        shape = eigen_problem.eigensolution.modeShapes[index]
        frequency = abs(shape.frequency)
        growth_rate = shape.damping
        damping_ratio = frequency > 0 ? -100 * growth_rate / hypot(growth_rate,frequency) : NaN
        (; index,eigenvalue=λ,frequency_rad_s=frequency,
            frequency_hz=frequency / (2π),growth_rate_per_s=growth_rate,
            damping_ratio_percent=damping_ratio,mode_shape=shape)
    end

    return (; analysis=:structural_fixed_aerodynamic_loads,eigen_problem,modes)
end

"""Write the fixed-load structural modal summary to CSV."""
function save_pazy_structural_modes(path::AbstractString,analysis)
    mkpath(dirname(path))
    open(path,"w") do io
        println(io,"analysis,mode,eigenvalue_real,eigenvalue_imag,frequency_hz,growth_rate_per_s,damping_ratio_percent")
        for mode in analysis.modes
            println(io,join((string(analysis.analysis),string(mode.index),
                string(real(mode.eigenvalue)),string(imag(mode.eigenvalue)),
                string(mode.frequency_hz),string(mode.growth_rate_per_s),
                string(mode.damping_ratio_percent)),','))
        end
    end
    return path
end

end