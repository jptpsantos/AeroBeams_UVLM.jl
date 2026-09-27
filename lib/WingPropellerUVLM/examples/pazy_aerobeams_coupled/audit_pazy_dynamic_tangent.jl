# Diagnostic only: test the discrete-time rotational chain rule against
# independently finite-differenced residuals. Production files are unchanged.
include(joinpath(@__DIR__, "replay_pazy_failure.jl"))

function add_time_chain_rule!(p)
    for e in p.model.elements
        scale = e.Δℓ / p.Δt / p.model.forceScaling
        dF = scale * AeroBeams.mul3(e.R_p1,e.R_p2,e.R_p3,e.R0*e.compStates.P)
        dM = scale * AeroBeams.mul3(e.R_p1,e.R_p2,e.R_p3,e.R0*e.compStates.H)
        p.jacobian[e.eqs_Fu1,e.DOF_p] .+= dF
        p.jacobian[e.eqs_Fu2,e.DOF_p] .+= dF
        p.jacobian[e.eqs_Fp1,e.DOF_p] .+= e.notHingedNode1Mat*dM
        p.jacobian[e.eqs_Fp2,e.DOF_p] .+= dM
    end
end

function tangent_audit(data)
    BLAS.set_num_threads(1)
    ctx = data.before
    p0 = deepcopy(ctx.structure)
    AeroBeams.get_equivalent_states_rates!(p0)
    offsets = 2 .* ctx.offsets .- ctx.previous_offsets
    loads = beam_loads(UVLM.dimensional_loads(data.aerodynamic.workspace.system), offsets, ctx.weights)
    apply_pazy_nodal_loads!(p0,loads,p0.timeNow)
    p0.skipJacobianUpdate = false
    AeroBeams.assemble_system_arrays!(p0)
    println("Failing physical time: ",p0.timeNow)
    println("Maximum element rotation-parameter norm: ",maximum(norm(e.states.p) for e in p0.model.elements))
    original = copy(p0.jacobian)
    add_time_chain_rule!(p0)
    corrected = copy(p0.jacobian)
    ids = unique(vcat([e.DOF_p for e in p0.model.elements]...))
    for seed in 1:3
        direction = zeros(length(p0.x))
        direction[ids] .= sin.(seed .* (1:length(ids)))
        direction ./= norm(direction)
        for h in (1e-5,1e-6,1e-7)
            plus,minus = deepcopy(p0),deepcopy(p0)
            AeroBeams.assemble_system_arrays!(plus,p0.x .+ h.*direction)
            AeroBeams.assemble_system_arrays!(minus,p0.x .- h.*direction)
            fd = (plus.residual-minus.residual)/(2h)
            println("Direction ",seed," h=",h," original error=",norm(fd-original*direction)/norm(fd),
                " corrected error=",norm(fd-corrected*direction)/norm(fd))
        end
    end
    # Reassemble at EVERY iteration, with identical initial state and loads.
    # Use actual residual convergence (not the solver's OR relative test).
    for (use_correction,backtrack) in ((false,false),(true,false),(true,true))
        p = deepcopy(p0)
        println("Newton with fresh tangent, chain-rule correction=",use_correction," backtracking=",backtrack)
        limit = backtrack ? 60 : 20
        for k in 0:limit
            p.skipJacobianUpdate = false
            AeroBeams.assemble_system_arrays!(p)
            r = norm(p.residual)
            println("  iteration=",k," reassembled residual=",r)
            r < 1e-6 && break
            k == limit && break
            use_correction && add_time_chain_rule!(p)
            increment = -(p.jacobian \ p.residual)
            if backtrack
                accepted = false
                for j in 0:20
                    fraction = 2.0^(-j)
                    trial = deepcopy(p)
                    trial.x .+= fraction .* increment
                    AeroBeams.assemble_system_arrays!(trial)
                    if norm(trial.residual) <= (1-1e-4*fraction)*r
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
                p.x .+= increment
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
