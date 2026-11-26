using AeroBeams, LinearAlgebra, DelimitedFiles
BLAS.set_num_threads(1)

function audit_case(io,sweep,pitch,mesh,rod,loss,airfoil,update,config,mass,rootCorr)
    Λ = deg2rad(sweep)
    θ = deg2rad(pitch + (sweep == 20 ? 0.4 : 0.0))
    _,L,c,spar,nodes = geometrical_properties_Pazy(; Λ,mesh)
    η = [0.0,-c*(1-spar) - 0.04,0.0]
    x0 = Float64[]
    for U in (0.001,10.0,20.0,30.0,40.0,50.0,60.0,70.0)
        sweep == 0 && U > 60 && break
        model,n = create_Pazy(; mesh,sweepAngle=Λ,pitchAngle=θ,airspeed=U,upright=true,
            beamWingtipRod=rod,withWingtipRod=true,sweptTipMassConfig=config,rootTorsionCorr=rootCorr,
            tipLossType=loss,hasTipCorrection=true,tipLossFunctionIsAirspeedDependent=update,
            aeroSolver=Indicial(),airfoil=deepcopy(airfoil),tipMasses=[mass],ηtipMasses=[η],
            tipMassesInertias=[zeros(3,3)])
        problem = create_SteadyProblem(; model,x0)
        solve!(problem)
        @assert problem.systemSolver.loadStepState.convergedFinalSolution
        x0 = copy(problem.systemArrays.x)
        wing = model.beams[1]
        tip = wing.elements[end].state.nodal
        gauge = argmin(abs.((nodes[1:end-1] .+ nodes[2:end])./2 .- 0.4))
        ε = wing.elements[gauge].state.midpoint
        strain = 1e6*(ε.γ[1] - ε.κ[2]*2.23e-3/2 + ε.κ[3]*50e-3/2)
        α = rad2deg(problem.loadstepHistory.aeroVariables[end][n].unsteadyKin.αₑ)
        dist = wing.elements[gauge].aero.loadsDistribution
        factor = tip_loss_function_Pazy(loss,θ,U,Λ)((nodes[gauge]+nodes[gauge+1])/2)
        @assert dist.Uref ≈ U
        @assert dist.ϖMid ≈ factor
        writedlm(io,permutedims(Any[sweep,pitch,U,mesh,Int(rod),loss,airfoil.name,Int(update),
            config,mass,Int(rootCorr),100tip.u2_b[3]/L,-100tip.u2[1]/L,twist_angle(tip.p2_b),
            strain,α,dist.ϖMid]),',')
        flush(io)
    end
end

open(joinpath(@__DIR__,"pazy_vlm_def_audit.csv"),"w") do io
    println(io,"sweep,pitch,speed,mesh,rod,loss,airfoil,update,config,mass,rootCorr,oop_local,oop_tunnel,twist,strain,aoa,factor")
    for pitch in (3,5,7), loss in ("VLM-def","Exponential")
        audit_case(io,0,pitch,:nominal,true,loss,FlatPlate,false,:baseline,0.01,true)
        audit_case(io,0,pitch,:refined,false,loss,NACA0018,false,:baseline,0.01,true)
        loss == "VLM-def" && audit_case(io,0,pitch,:refined,true,loss,FlatPlate,false,:baseline,0.01,true)
    end
    for sweep in (10,20), pitch in (3,5,7,10), mesh in (:nominal,:refined)
        audit_case(io,sweep,pitch,mesh,false,"VLM-def",FlatPlate,false,:TE,0.0,true)
    end
end
