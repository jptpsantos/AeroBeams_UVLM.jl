# Exercise the saver directly with a result created before the local
# perturbation record was added, then with the new record and local-only mode.
run_pazy_wing_uvlm_static = nothing
include(joinpath(@__DIR__,"..","lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled","PazyWingUVLMStaticSweep.jl"))
z = [0.0]
section = (; cn=z,cm_quarter_chord=z,cm_spar=z,eta=z,s_m=z,y_reference_m=z,
    positions_UVLM=zeros(3,1),alpha_local=z)
oldCase = (; section,airspeed=50.0,angle_of_attack=0.0,sweepAngle=0.0,alpha_step=deg2rad(0.25),
    cn_alpha_equilibrium=z,cn_alpha_fixed_shape=[4.0],central_converged=true,
    equilibrium_slope_converged=false,fixed_shape_slope_converged=true,
    tip_out_of_plane=0.0,tip_twist_degrees=0.0,cn_alpha_local=[NaN],cn_alpha_local_valid=[false],
    cm_quarter_chord_alpha_equilibrium=z,cm_spar_alpha_equilibrium=z,
    cm_quarter_chord_alpha_fixed_shape=z,cm_spar_alpha_fixed_shape=z,
    diagnostics=[(; solve=:central,converged=true,reason="Mock central solve")])
result = (; cases=[oldCase],slope_modes=(:fixed_shape,),static_options=(;),
    local_alpha_tolerance=1e-8,maximum_aerodynamic_chords=80.0,
    minimum_aerodynamic_chords=20.0,aerodynamic_window_chords=4.0)
mktempdir(@__DIR__) do directory
    path = joinpath(directory,"aoa0deg","spanwise_cn_alpha_perturbations.csv")
    save_pazy_static_sweep(result; output_directory=directory,save_plots=false)
    table = readdlm(path,',',Any)
    @assert size(table) == (2,13) && all(isnan,table[2,10:12]) && string(table[2,13]) == "false"
    fixedShape = (; delta_cn=[0.02],delta_alpha_rad=[0.005],cn_alpha=[4.0],valid=[true])
    current = merge(oldCase,(; cn_alpha_local_perturbations=(; fixed_shape=fixedShape)))
    save_pazy_static_sweep(merge(result,(; cases=[current])); output_directory=directory,save_plots=false)
    table = readdlm(path,',',Any)
    @assert table[2,12] == 4.0 && string(table[2,13]) == "true"
    save_pazy_static_sweep(merge(result,(; slope_modes=(:local_ratio,)));
        output_directory=directory,save_plots=false)
    @assert size(readdlm(path,',',Any)) == (1,13)
end
println("Legacy results, current perturbation results, and local-only overwrite verified")
