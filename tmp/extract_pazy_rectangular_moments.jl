using LinearAlgebra
BLAS.set_num_threads(1)
include(joinpath(@__DIR__,"..","lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled","PazyWingUVLMStaticSweep.jl"))

# Keep a matched archived-mesh zero-incidence case and use a smaller lattice
# for the two deflected cases.
archivedMesh = "--archived-mesh" in ARGS
nc,ns = archivedMesh ? (5,15) : (3,8)
options = (; mesh=:nominal,density=1.225,chordwise_panels=nc,spanwise_panels=ns,
    symmetric_wing=true,time_step_chords=1 / nc,maximum_wake_rows=10nc,core_radius=1e-3,
    wake_shedding_fraction=0.1,aerodynamic_interaction=true,
    aerodynamic_force_tolerance=1e-5,aerodynamic_circulation_tolerance=1e-5,
    maximum_static_iterations=100,static_relaxation=archivedMesh ? 0.3 : 0.55,geometry_tolerance=1e-4,
    load_tolerance=1e-4,consecutive_equilibrium_iterations=2,
    newton_maximum_iterations=50,newton_absolute_tolerance=1e-8,
    newton_relative_tolerance=1e-8,newton_display_iterations=false,progress_frequency=1000)
outputDirectory = joinpath(@__DIR__,"..","lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled","output",
    "static_spanwise_cn_cm","rectangular_diagnostic")
if !archivedMesh
    outputDirectory = joinpath(outputDirectory,"reduced_uvlm_mesh")
end
for (speed,incidence) in (archivedMesh ? ((60.0,0.0),) : ((50.0,3.0),(120.0,1.0)))
    result = run_pazy_wing_uvlm_static_sweep(;
        airspeeds=[speed],angles_of_attack=deg2rad.([incidence]),
        static_options=options,slope_modes=(:local_ratio,:fixed_shape),output_directory=outputDirectory)
    case = only(result.cases)
    @assert case.central_converged && case.fixed_shape_slope_converged
    j = 2
    println("RESULT U=$speed alpha=$incidence eta=$(case.section.eta[j]) " *
        "cn=$(case.section.cn[j]) secant=$(case.cn_alpha_local[j]) " *
        "fixed_shape=$(case.cn_alpha_fixed_shape[j]) " *
        "cm_quarter=$(case.section.cm_quarter_chord[j]) cm_spar=$(case.section.cm_spar[j]) " *
        "cm_spar_slope=$(case.cm_spar_alpha_fixed_shape[j])")
    flush(stdout)
end
