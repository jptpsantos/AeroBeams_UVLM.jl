using LinearAlgebra
BLAS.set_num_threads(1)
const exampleDir = joinpath(@__DIR__,"..","lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled")
include(joinpath(exampleDir,"test","test_static_spanwise_cn_alpha.jl"))

directory = joinpath(exampleDir,"output","static_spanwise_cn_cm","perturbation_verification")
result = run_pazy_wing_uvlm_static_sweep(; airspeeds=[50.0],angles_of_attack=deg2rad.([3.0]),
    alpha_step=deg2rad(0.25),slope_modes=(:local_ratio,:fixed_shape,:equilibrium),
    static_options=options,maximum_aerodynamic_chords=80.0,output_directory=directory,save_plots=false)
case = only(result.cases)
@assert all(case.cn_alpha_local_perturbations.equilibrium.valid)
@assert all(case.cn_alpha_local_perturbations.fixed_shape.valid)
j = 2
println("EXPORT ",(; eta=case.section.eta[j],ratio=case.cn_alpha_local[j],
    equilibrium_freestream=case.cn_alpha_equilibrium[j],fixed_shape_freestream=case.cn_alpha_fixed_shape[j],
    equilibrium_local=case.cn_alpha_local_perturbations.equilibrium.cn_alpha[j],
    fixed_shape_local=case.cn_alpha_local_perturbations.fixed_shape.cn_alpha[j],
    equilibrium_delta_local_deg=rad2deg(case.cn_alpha_local_perturbations.equilibrium.delta_alpha_rad[j]),
    fixed_shape_delta_local_deg=rad2deg(case.cn_alpha_local_perturbations.fixed_shape.delta_alpha_rad[j])))
