# Replay the SAME failing physical step from audit_pazy_failure.jl's checkpoint.
# No time advance, load clipping, or changes to the production solvers.
using LinearAlgebra, Serialization
include(joinpath(@__DIR__, "PazyWingUVLMCoupling.jl"))

function replay_failure_data(data)
    ctx = data.before
    @assert ctx.i == data.context.i "Checkpoint does not precede the failing step"
    p0 = deepcopy(ctx.structure)
    AeroBeams.get_equivalent_states_rates!(p0)
    offsets = 2 .* ctx.offsets .- ctx.previous_offsets
    loads = beam_loads(UVLM.dimensional_loads(data.aerodynamic.workspace.system), offsets, ctx.weights)
    apply_pazy_nodal_loads!(p0, loads, p0.timeNow)
    println("Replaying t=", p0.timeNow, ", dt=", p0.Δt)
    println("Max nodal force=", maximum(norm, eachcol(loads[1:3,:])),
        ", max nodal moment=", maximum(norm, eachcol(loads[4:6,:])))

    for (name, fresh, loadscale) in (("original", false, 1.0), ("fresh Jacobian", true, 1.0),
                                   ("zero applied loads", true, 0.0))
        p = deepcopy(p0)
        p.systemSolver.alwaysUpdateJacobian = fresh
        p.systemSolver.displayStatus = false
        apply_pazy_nodal_loads!(p, loadscale .* loads, p.timeNow)
        p.skipJacobianUpdate = false
        AeroBeams.assemble_system_arrays!(p)
        first_residual = norm(p.residual)
        AeroBeams.solve_time_step!(p)
        AeroBeams.assemble_system_arrays!(p)
        println(name, ": converged=", p.systemSolver.convergedFinalSolution,
            ", initial residual=", first_residual, ", final reassembled residual=", norm(p.residual))
    end

    p0.skipJacobianUpdate = false
    AeroBeams.assemble_system_arrays!(p0)
    J = copy(p0.jacobian)
    x = copy(p0.x)
    println("Jacobian condition estimate (2-norm)=", cond(Matrix(J)))
    for block in (:DOF_u, :DOF_p, :DOF_V, :DOF_Ω, :DOF_F, :DOF_M)
        ids = unique(vcat([getproperty(e, block) for e in p0.model.elements]...))
        direction = zeros(length(x))
        direction[ids] .= sin.(1:length(ids))
        direction ./= norm(direction)
        for h in (1e-5, 1e-6, 1e-7)
            plus, minus = deepcopy(p0), deepcopy(p0)
            AeroBeams.assemble_system_arrays!(plus, x .+ h .* direction)
            AeroBeams.assemble_system_arrays!(minus, x .- h .* direction)
            fd = (plus.residual-minus.residual)/(2h)
            println("Jacobian ", block, " h=", h, " relative directional error=",
                norm(fd-J*direction)/max(norm(fd), eps()))
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    checkpoint = only(ARGS)
    data = deserialize(checkpoint)
    # Deserializing user-defined motion functions can introduce new methods.
    # Enter their world before calling the structural assembly routines.
    path = joinpath(dirname(checkpoint), "step_replay.txt")
    open(path, "w") do io
        redirect_stdout(io) do
            Base.invokelatest(replay_failure_data, data)
        end
    end
    print(read(path, String))
end
