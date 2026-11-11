using Test,LinearAlgebra
BLAS.set_num_threads(1)
ENV["GKSwstype"] = "100"

@testset "Additive AeroBeams wing-propeller coupling" begin
    for relative in ("../test_reduced_rotor.jl","../test_tip_rotor.jl",
        "test_aerodynamic_interface.jl","test_coupled_solver.jl",
        "check_pylon.jl","check_mean_equilibrium.jl")
        # Isolate standalone verification names while sharing loaded packages.
        scope = Module(gensym(:ChangVerification))
        Core.eval(scope,:(include(path::AbstractString) = Base.include($scope,path)))
        Base.include(scope,joinpath(@__DIR__,relative))
    end
end
