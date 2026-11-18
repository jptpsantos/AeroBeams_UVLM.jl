using LinearAlgebra
BLAS.set_num_threads(1)
const exampleDir = joinpath(@__DIR__,"..","lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled")
include(joinpath(exampleDir,"test","test_static_spanwise_cn_alpha.jl"))

pilotOptions = (; options...,mesh=:nominal,chordwise_panels=3,spanwise_panels=8,
    time_step_chords=1/3,maximum_wake_rows=30,static_relaxation=0.3,
    load_tolerance=1e-4,load_residual_force_floor=1e-9)
directory = joinpath(exampleDir,"output","static_spanwise_cn_cm","near_zero_verification")
result = run_pazy_wing_uvlm_static_sweep(; airspeeds=[60.0],angles_of_attack=deg2rad.([1e-6]),
    local_alpha_tolerance=1e-10,static_options=pilotOptions,output_directory=directory)
case = only(result.cases)
@assert case.central_converged && all(case.cn_alpha_local_valid)
println("EXPORT ",(; folder=pazy_static_angle_folder(case.angle_of_attack),
    eta=case.section.eta[2],alpha_deg=rad2deg(case.section.alpha_local[2]),
    cn=case.section.cn[2],cn_alpha_local=case.cn_alpha_local[2],
    cm_quarter_chord=case.section.cm_quarter_chord[2],cm_spar=case.section.cm_spar[2]))
