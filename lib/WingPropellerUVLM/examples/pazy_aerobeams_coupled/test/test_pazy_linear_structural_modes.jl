using Test
include(joinpath(@__DIR__,"..","PazyWingUVLMStatic.jl"))
include(joinpath(@__DIR__,"..","PazyLinearStructuralModes.jl"))

@testset "Pazy fixed-load structural modes" begin
    @test_throws ArgumentError PazyLinearStructuralModes.solve_pazy_structural_modes((;))
    @test_throws ArgumentError PazyLinearStructuralModes.solve_pazy_structural_modes(
        (;converged=false,structural_problem=nothing))

    pazy = create_pazy_uvlm_structure(;mesh=:nominal,withWingtipRod=false)
    solver = AeroBeams.create_NLsolver(method=:newton,maxIter=30,ftol=1e-8,
        displayStatus=false,trackLoadSteps=false)
    structural_problem = AeroBeams.create_SteadyProblem(
        model=pazy.model,systemSolver=solver)
    AeroBeams.solve!(structural_problem)
    @test structural_problem.systemSolver.loadStepState.convergedFinalSolution
    reference = (;converged=true,structural_problem)
    reference_state = copy(structural_problem.systemArrays.x)

    analysis = PazyLinearStructuralModes.solve_pazy_structural_modes(
        reference;number_of_modes=4,eigenmethod=:fulldecomp)
    @test analysis.analysis == :structural_fixed_aerodynamic_loads
    @test !isempty(analysis.modes)
    @test all(mode -> isfinite(mode.frequency_hz) && mode.frequency_hz >= 0,
        analysis.modes)
    @test structural_problem.systemArrays.x == reference_state

    mktempdir() do directory
        path = PazyLinearStructuralModes.save_pazy_structural_modes(
            joinpath(directory,"modes.csv"),analysis)
        @test isfile(path)
        @test occursin("structural_fixed_aerodynamic_loads",read(path,String))
    end
end