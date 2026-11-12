# Diagnostic only: test the discrete-time rotational chain rule against
# independently finite-differenced residuals. Production files are unchanged.
include(joinpath(@__DIR__,"replay_pazy_failure.jl"))

function add_time_chain_rule!(p)
    for e in p.model.assemblyData.elements
        scale = e.geometry.Δℓ / p.time.state.Δt / p.model.assemblyData.forceScaling
        dF = scale * AeroBeams.mul3(e.stateCache.rotationMatrices.R_p1,e.stateCache.rotationMatrices.R_p2,e.stateCache.rotationMatrices.R_p3,e.geometry.R0 * e.state.midpoint.P)
        dM = scale * AeroBeams.mul3(e.stateCache.rotationMatrices.R_p1,e.stateCache.rotationMatrices.R_p2,e.stateCache.rotationMatrices.R_p3,e.geometry.R0 * e.state.midpoint.H)
        p.systemArrays.jacobian[e.assemblyIndices.eqs_Fu1,e.assemblyIndices.DOF_p] .+= dF
        p.systemArrays.jacobian[e.assemblyIndices.eqs_Fu2,e.assemblyIndices.DOF_p] .+= dF
        p.systemArrays.jacobian[e.assemblyIndices.eqs_Fp1,e.assemblyIndices.DOF_p] .+= e.geometry.notHingedNode1Mat * dM
        p.systemArrays.jacobian[e.assemblyIndices.eqs_Fp2,e.assemblyIndices.DOF_p] .+= dM
    end
end

function tangent_audit(data)
    BLAS.set_num_threads(1)
    ctx = data.before
    p0 = deepcopy(ctx.structure)
    AeroBeams.update_equivalent_states_rates!(p0)
    offsets = 2 .* ctx.offsets .- ctx.previous_offsets
    loads = beam_loads(UVLM.dimensional_loads(data.aerodynamic.workspace.system),offsets,ctx.weights)
    apply_pazy_nodal_loads!(p0,loads,p0.time.state.timeNow)
    p0.systemSolver.state.skipJacobianUpdate = false
    AeroBeams.fj!(p0)
    println("Failing physical time: ",p0.time.state.timeNow)
    println("Maximum element rotation-parameter norm: ",maximum(norm(e.state.midpoint.p) for e in p0.model.assemblyData.elements))
    original = copy(p0.systemArrays.jacobian)
    add_time_chain_rule!(p0)
    corrected = copy(p0.systemArrays.jacobian)
    ids = unique(vcat([e.assemblyIndices.DOF_p for e in p0.model.assemblyData.elements]...))
    for seed = 1:3
        direction = zeros(length(p0.systemArrays.x))
        direction[ids] .= sin.(seed .* (1:length(ids)))
        direction ./= norm(direction)
        for h in (1e-5,1e-6,1e-7)
            plus,minus = deepcopy(p0),deepcopy(p0)
            plus.systemArrays.x .= p0.systemArrays.x .+ h .* direction
            minus.systemArrays.x .= p0.systemArrays.x .- h .* direction
            AeroBeams.fj!(plus)
            AeroBeams.fj!(minus)
            fd = (plus.systemArrays.residual - minus.systemArrays.residual) / (2h)
            println("Direction ",seed," h=",h," original error=",norm(fd - original * direction) / norm(fd),
                " corrected error=",norm(fd - corrected * direction) / norm(fd))
        end
    end
    # Reassemble at EVERY iteration, with identical initial state and loads.
    # Use actual residual convergence (not the solver's OR relative test).
    for (useCorrection,backtrack) in ((false,false),(true,false),(true,true))
        p = deepcopy(p0)
        println("Newton with fresh tangent, chain-rule correction=",useCorrection," backtracking=",backtrack)
        limit = backtrack ? 60 : 20
        for k = 0:limit
            p.systemSolver.state.skipJacobianUpdate = false
            AeroBeams.fj!(p)
            r = norm(p.systemArrays.residual)
            println("  iteration=",k," reassembled residual=",r)
            r < 1e-6 && break
            k == limit && break
            useCorrection && add_time_chain_rule!(p)
            increment = -(p.systemArrays.jacobian \ p.systemArrays.residual)
            if backtrack
                accepted = false
                for j = 0:20
                    fraction = 2.0^(-j)
                    trial = deepcopy(p)
                    trial.systemArrays.x .+= fraction .* increment
                    AeroBeams.fj!(trial)
                    if norm(trial.systemArrays.residual) <= (1 - 1e-4 * fraction) * r
                        println("    accepted Newton fraction=",fraction)
                        p = trial
                        accepted = true
                        break
                    end
                end
                if !accepted
                    println("    backtracking failed to reduce residual")
                    break
                end
            else
                p.systemArrays.x .+= increment
            end
        end
    end
end

data = deserialize(only(ARGS))
path = joinpath(dirname(only(ARGS)),"dynamic_tangent_audit.txt")
open(path,"w") do io
    redirect_stdout(io) do
        Base.invokelatest(tangent_audit,data)
    end
end
print(read(path,String))
