using Test
using DelimitedFiles

# Exercise the actual output functions without rerunning the unchanged FSI solver.
sourcePath = joinpath(@__DIR__,"..","lib","WingPropellerUVLM","examples",
    "pazy_aerobeams_coupled","PazyWingUVLMStaticSweep.jl")
definitions = Meta.parseall(read(sourcePath,String))
for definition in definitions.args
    if definition isa Expr && definition.head == :function
        name = definition.args[1].args[1]
        if name in (:pazy_static_angle_folder,:save_pazy_static_sweep,:save_pazy_static_sweep_cases)
            Core.eval(@__MODULE__,definition)
        end
    end
end
if "--plots" in ARGS
    import Plots
    const AeroBeams = (; Plots)
end

function output_case(angle,U;sweep=0.0)
    section = (; eta=[0.25,0.75],s_m=[0.15,0.45],y_reference_m=[0.15,0.45],
        positions_UVLM=zeros(3,2),cn=[0.2,0.1],alpha_local=fill(deg2rad(angle),2))
    diagnostic = (; solve=:central,converged=true,reason="Static equilibrium reached",
        history=[U 0.0],history_columns=["U","residual"])
    return (; airspeed=U,angle_of_attack=deg2rad(angle),sweepAngle=deg2rad(sweep),alpha_step=NaN,
        section,cn_alpha_local=section.cn ./ section.alpha_local,cn_alpha_local_valid=trues(2),
        cn_alpha_equilibrium=fill(NaN,2),cn_alpha_fixed_shape=fill(NaN,2),central_converged=true,
        equilibrium_slope_converged=false,fixed_shape_slope_converged=false,
        tip_out_of_plane=0.0,tip_twist_degrees=0.0,diagnostics=[diagnostic])
end

directory = mktempdir(@__DIR__; cleanup=false)
result = (; cases=[output_case(3.0,10.0),output_case(-2.5,20.0),output_case(3.0,30.0;sweep=20.0)],
    slope_modes=(:local_ratio,),static_options=(;),local_alpha_tolerance=1e-8,
    maximum_aerodynamic_chords=80.0,minimum_aerodynamic_chords=25.0,aerodynamic_window_chords=4.0)

@testset "Output folders by prescribed incidence" begin
    @test pazy_static_angle_folder(deg2rad(3.0)) == "aoa3deg"
    @test pazy_static_angle_folder(deg2rad(2.5)) == "aoa2p5deg"
    @test pazy_static_angle_folder(deg2rad(-2.5)) == "aoam2p5deg"
    @test pazy_static_angle_folder(0.0) == "aoa0deg"
    @test pazy_static_angle_folder(deg2rad(3.0000000001)) != pazy_static_angle_folder(deg2rad(3.0))
    # Incremental checkpoint first, followed by a multi-angle result.
    save_pazy_static_sweep(merge(result,(;cases=result.cases[1:1])); output_directory=directory,save_plots=false)
    @test save_pazy_static_sweep(result; output_directory=directory,save_plots="--plots" in ARGS) == directory
    @test !isfile(joinpath(directory,"spanwise_cn_alpha.csv"))
    @test Set(readdir(directory)) == Set(["aoa3deg","aoam2p5deg"])
    a = readdlm(joinpath(directory,"aoa3deg","spanwise_cn_alpha.csv"),',',Any)
    b = readdlm(joinpath(directory,"aoam2p5deg","spanwise_cn_alpha.csv"),',',Any)
    @test size(a) == (5,23)
    @test size(b) == (3,23)
    @test all(x -> isapprox(x,3.0),a[2:end,3])
    @test all(x -> isapprox(x,-2.5),b[2:end,3])
    @test a[2:end,1] == [1,1,2,2]
    @test Set(a[2:end,2]) == Set([10.0,30.0])
    @test Set(b[2:end,2]) == Set([20.0])
    @test isfile(joinpath(directory,"aoa3deg","case_2_central_convergence.csv"))
    @test readdlm(joinpath(directory,"aoa3deg","case_2_central_convergence.csv"),',',Any)[2,1] == 30.0
    @test size(readdlm(joinpath(directory,"aoa3deg","status.csv"),',',Any),1) == 3
    @test size(readdlm(joinpath(directory,"aoam2p5deg","status.csv"),',',Any),1) == 2
    @test isfile(joinpath(directory,"aoam2p5deg","settings.txt"))
    original = read(joinpath(directory,"aoa3deg","spanwise_cn_alpha.csv"),String)
    save_pazy_static_sweep(merge(result,(;cases=result.cases[2:2])); output_directory=directory,save_plots=false)
    @test read(joinpath(directory,"aoa3deg","spanwise_cn_alpha.csv"),String) == original
    if "--plots" in ARGS
        @test count(file -> endswith(file,".png"),readdir(joinpath(directory,"aoa3deg"))) == 2
        @test count(file -> endswith(file,".png"),readdir(joinpath(directory,"aoam2p5deg"))) == 1
    end
end
println("Angle-folder output verification: ",directory)
