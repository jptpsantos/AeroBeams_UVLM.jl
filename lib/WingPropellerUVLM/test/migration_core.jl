using Test, WingPropellerUVLM, LinearAlgebra
import WingPropellerUVLM as U

@testset "Model/problem migration and accepted wake history" begin
    grid,_=wing_to_grid([0.,0.],[0.,1.],[0.,0.],ones(2),zeros(2),zeros(2),2,2)
    model=create_UVLMModel(;surfaces=[grid],reference=Reference(1.,1.,1.,zeros(3),10.))
    op=create_OperatingPoint(;airspeed=10.,angleOfAttack=deg2rad(3.))
    solver=create_UVLMSolver(;coreRadius=.001)
    makeproblem()=create_UVLMDynamicProblem(;model,operatingPoint=op,aeroSolver=solver,
        timeVector=[0.,.01,.02,.03],maximumWakeRows=2,saveGeometry=true)
    p=makeproblem()
    @test p.state.system !== p.workspace.system
    @test p.state.system.Γ !== p.workspace.system.Γ
    @test p.state.timeNow==0.
    @test_throws ArgumentError commit_time_step!(p)
    @test_throws ArgumentError evaluate_trial!(p)
    for step in 1:3
        # Independent legacy one-call step is the numerical reference, including
        # mature-wake convection and capacity truncation on the third step.
        reference=deepcopy(p.state.system)
        U.copy_surfaces_to_previous!(reference,1)
        propagate_system!(reference,U.freestream(op),.01;
            additional_velocity=nothing,repeated_points=U.repeated_trailing_edge_points(reference.surfaces),
            nwake=copy(reference.nwake),eta=.1,calculate_influence_matrix=true,
            near_field_analysis=true,derivatives=false,advance_wake=true)
        begin_time_step!(p)
        evaluate_trial!(p)
        firstGamma=copy(p.workspace.system.Γ)
        displaced=deepcopy(model.grids); displaced[1][3,:,:] .+= .002
        evaluate_trial!(p;surfaces=displaced)
        evaluate_trial!(p)
        @test p.workspace.system.Γ==firstGamma
        @test p.state.stepIndex==step
        @test length(p.results.savedTimeVector)==step-1
        @test p.workspace.system.Γ ≈ reference.Γ
        @test p.workspace.system.dΓdt ≈ reference.dΓdt
        @test p.workspace.system.chord_seg_forces ≈ reference.chord_seg_forces
        commit_time_step!(p)
        @test all(isapprox(getfield(a,k),getfield(b,k)) for
            (a,b) in zip(p.state.system.wakes[1][1:min(step,2),:],reference.wakes[1][1:min(step,2),:]) for
            k in fieldnames(typeof(a)))
        @test p.state.system.nwake==[min(step,2)]
        @test_throws ArgumentError commit_time_step!(p)
    end
    @test p.results.completed
    @test p.results.savedTimeVector==[.01,.02,.03]
    @test all(isfinite,reduce(vcat,p.results.forceOverTime))
    saved=copy(p.results.circulationOverTime[end])
    p.state.system.Γ .= 0
    @test p.results.circulationOverTime[end]==saved
    @test model.grids[1]==grid
    @test solve!(makeproblem()).results.completed
    bad=makeproblem(); begin_time_step!(bad)
    invalid=deepcopy(model.grids); invalid[1][1]=NaN
    @test_throws ArgumentError evaluate_trial!(bad;surfaces=invalid)
    @test_throws ArgumentError commit_time_step!(bad)
    rollback_time_step!(bad)
    @test bad.state.timeNow==0 && isempty(bad.results.savedTimeVector)
    @test solve!(bad).results.completed
    steady=create_UVLMSteadyProblem(;model,operatingPoint=op,aeroSolver=solver)
    solve!(steady)
    direct=steady_analysis([grid],Reference(1.,1.,1.,zeros(3),10.,op.density),U.freestream(op);
        fcore=(c,ds)->.001,derivatives=false)
    @test steady.system.Γ ≈ direct.Γ
    scaled=create_OperatingPoint(;airspeed=85.,referenceAirspeed=65.,referenceRotationRPM=1200.,
        rotorSpeedPolicy=:constantAdvanceRatio)
    @test scaled.angularSpeed ≈ 2pi/60*1200*85/65
    @test_throws ArgumentError create_OperatingPoint(;airspeed=-1.)
    @test_throws ArgumentError create_OperatingPoint(;airspeed=85.,rotationRPM=12.,referenceRotationRPM=12.,rotorSpeedPolicy=:constantAdvanceRatio)
    @test_throws ArgumentError create_UVLMDynamicProblem(;model,operatingPoint=op,timeVector=[0.,0.])
    @test_throws ArgumentError create_UVLMSolver(;coreRadius=0.)
    @test create_PartitionedCoupling(;maximumIterations=3).maximum_iterations==3
end

