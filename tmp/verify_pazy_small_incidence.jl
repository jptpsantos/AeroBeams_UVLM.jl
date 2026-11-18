using LinearAlgebra
BLAS.set_num_threads(1)
include(joinpath(@__DIR__,"..","lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled","PazyWingUVLMStaticSweep.jl"))

options = (; mesh=:nominal,density=1.225,chordwise_panels=2,spanwise_panels=3,
    symmetric_wing=true,time_step_chords=0.5,maximum_wake_rows=6,core_radius=1e-3,
    wake_shedding_fraction=0.1,aerodynamic_interaction=true,
    aerodynamic_force_tolerance=1e-5,aerodynamic_circulation_tolerance=1e-5,
    maximum_static_iterations=150,static_relaxation=0.55,geometry_tolerance=1e-4,
    load_tolerance=1e-4,consecutive_equilibrium_iterations=2,
    newton_maximum_iterations=50,newton_absolute_tolerance=1e-8,
    newton_relative_tolerance=1e-8,newton_display_iterations=false,progress_frequency=1000)

if "--runner-tolerances" in ARGS
    for speed in [100.0,120.0]
        result = run_pazy_wing_uvlm_static_sweep(;
            airspeeds=[speed],angles_of_attack=[deg2rad(1e-6)],
            static_options=(; options...,static_relaxation=0.3,maximum_static_iterations=100,
                load_residual_force_floor=1e-9),local_alpha_tolerance=1e-10)
        case = only(result.cases)
        println("RUNNER TOLERANCES ",(; speed,case.central_converged,case.tip_out_of_plane,
            ratio=case.cn_alpha_local,history=only(case.diagnostics).history[end,:]))
        flush(stdout)
    end
    exit()
end

for speed in [60.0,100.0],(incidence,tolerance) in [(1e-6,1e-4),(1e-6,1e-12),(1e-4,1e-10)]
    result = run_pazy_wing_uvlm_static_sweep(;
        airspeeds=[speed],angles_of_attack=deg2rad.([incidence]),
        static_options=(; options...,load_tolerance=tolerance),local_alpha_tolerance=1e-10)
    case = only(result.cases)
    println("RESULT ",(; speed,incidence,tolerance,case.central_converged,
        alpha=rad2deg.(case.section.alpha_local),ratio=case.cn_alpha_local,
        cm=case.section.cm_spar,tip=case.tip_out_of_plane))
    flush(stdout)
end
