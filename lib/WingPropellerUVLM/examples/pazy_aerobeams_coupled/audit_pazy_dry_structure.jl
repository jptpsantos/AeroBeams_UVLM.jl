# Isolate AeroBeams: same Pazy beam, time step and Newton settings, no UVLM.
# A short smooth force/moment pulse is followed by free vibration.
using LinearAlgebra, Dates
include(joinpath(@__DIR__, "PazyWingUVLMCoupling.jl"))

function energy(model)
    sum(model.elements) do e
        s, c = e.states, e.compStates
        e.Δℓ/2 * (dot(s.F,c.γ) + dot(s.M,c.κ) + dot(s.V,c.P) + dot(s.Ω,c.H))
    end
end

function main()
    BLAS.set_num_threads(1)
    ne, span, chord, spar = AeroBeams.geometrical_properties_Pazy()
    beam = AeroBeams.create_Beam(name="Dry Pazy audit", length=span, nElements=ne,
        normalizedNodalPositions=AeroBeams.nodal_positions_Pazy(),
        S=AeroBeams.stiffness_matrices_Pazy(withSkin=true, sweepStructuralCorrections=false, GAy=1e16, GAz=1e16, Λ=0.),
        I=AeroBeams.inertia_matrices_Pazy(withSkin=true),
        rotationParametrization="E231", p0=[-pi/2;0.;0.], aeroSurface=nothing)
    root = AeroBeams.create_BC(beam=beam, node=1,
        types=["u1A","u2A","u3A","p1A","p2A","p3A"], values=zeros(6))
    bcs = [root]
    for n in 1:ne+1
        push!(bcs, AeroBeams.create_BC(beam=beam,node=n,
            types=["F1A","F2A","F3A","M1A","M2A","M3A"],values=zeros(6)))
    end
    model = AeroBeams.create_Model(beams=[beam],BCs=bcs,gravityVector=zeros(3),v_A=t->zeros(3))
    dt = 0.25*chord/50
    times = collect(0:round(Int,3/dt)) .* dt
    nr = AeroBeams.create_NewtonRaphson(maximumIterations=20,absoluteTolerance=1e-6,
        relativeTolerance=1e-6,alwaysUpdateJacobian=false,minConvRateJacUpdate=1.2,displayStatus=false)
    p = AeroBeams.create_DynamicProblem(model=model,timeVector=times,systemSolver=nr,displayProgress=false)
    AeroBeams.precompute_distributed_loads!(p)
    AeroBeams.solve_initial_dynamic!(p)
    folder = joinpath(@__DIR__, "output", "dry_audit_" * Dates.format(now(), "yyyymmdd_HHMMSS"))
    mkpath(folder)
    energies = Float64[]
    loads = zeros(6,ne+1)
    open(joinpath(folder,"history.tsv"),"w") do io
        println(io,"time\tenergy_J\ttip_bending_m\ttip_twist_deg\tmax_omega_rad_s")
        for i in 2:length(times)
            t = times[i]
            AeroBeams.update_time_variables!(p,i)
            AeroBeams.update_basis_A_orientation!(p)
            AeroBeams.get_equivalent_states_rates!(p)
            pulse = 0 < t < 0.02 ? sinpi(t/0.02)^2 : 0.
            loads[1,end] = -pulse
            loads[6,end] = 0.03*pulse
            apply_pazy_nodal_loads!(p,loads,t)
            AeroBeams.solve_time_step!(p)
            p.systemSolver.convergedFinalSolution || error("Dry beam failed at t=$t")
            e = energy(p.model)
            t > 0.03 && push!(energies,e)
            println(io,join((t,e,-p.model.elements[end].nodalStates.u_n2[1],
                wingtip_twist_degrees(p.model),maximum(norm(e.states.Ω) for e in p.model.elements)), '\t'))
            if i % 1000 == 0
                flush(io)
                println("Dry audit t=",t," E=",e)
            end
        end
    end
    report = "Completed $(length(times)-1) dry steps to $(last(times)) s.\n" *
        "Free-vibration energy range: $(extrema(energies)) J\n" *
        "Final / initial free-vibration energy: $(last(energies)/first(energies))\n"
    open(joinpath(folder,"summary.txt"),"w") do io
        print(io,report)
    end
    println(report,"Saved ",folder)
end
main()
