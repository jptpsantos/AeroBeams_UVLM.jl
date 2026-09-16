# Run with --project=lib/WingPropellerUVLM. Including defines the module only.
module PrescribedWingExample
import WingPropellerUVLM as UVLM

function create_problem()
    grid,_=UVLM.wing_to_grid([0.,0.],[0.,1.],[0.,0.],ones(2),zeros(2),zeros(2),4,2)
    model=UVLM.create_UVLMModel(;surfaces=[grid],
        reference=UVLM.Reference(1.,1.,1.,zeros(3),10.),name="Prescribed heaving wing")
    point=UVLM.create_OperatingPoint(;airspeed=10.,density=1.225,angleOfAttack=deg2rad(3.))
    solver=UVLM.create_UVLMSolver(;coreRadius=.001)
    function motion(model,t)
        grids=deepcopy(model.grids)
        grids[1][3,:,:] .+= .002*sin(2pi*2t)
        return grids
    end
    return UVLM.create_UVLMDynamicProblem(;model,operatingPoint=point,aeroSolver=solver,
        timeVector=collect(0.:.005:.05),maximumWakeRows=10,surfaceMotion=motion,
        trackingTimeSteps=true,trackingFrequency=1,saveGeometry=false)
end
main()=UVLM.solve!(create_problem())
end
if abspath(PROGRAM_FILE)==@__FILE__
    problem=PrescribedWingExample.main()
    println("Completed: ",problem.results.completed)
    println("Last force [N]: ",last(problem.results.forceOverTime))
end
