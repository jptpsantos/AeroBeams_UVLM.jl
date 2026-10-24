# Diagnostic only. Same equilibrium and discretization for every case.
using LinearAlgebra, DelimitedFiles, Statistics, Printf, Serialization
BLAS.set_num_threads(1)
const BASE = normpath(joinpath(@__DIR__, "..", "lib", "WingPropellerUVLM", "examples", "pazy_aerobeams_coupled"))
include(joinpath(BASE,"PazyWingUVLMStatic.jl"))
const OUT = joinpath(BASE,"output","LCO_small_signal_audit")
mkpath(OUT)
const DT = 20 / ceil(Int,20/(0.25*0.0989/52))
const COMMON = (; airspeed=45.0,density=1.225,angle_of_attack=deg2rad(5.),
    chordwise_panels=5,spanwise_panels=16,symmetric_wing=true,
    maximum_wake_rows=18,time_step_chords=DT*45/0.0989,core_radius=1e-3,
    wake_shedding_fraction=0.1,aerodynamic_interaction=true,progress_frequency=500)
const REFINE_ONLY = "--refine" in ARGS
equilibrium = REFINE_ONLY ? deserialize(joinpath(OUT,"equilibrium.jls")) : run_pazy_wing_uvlm_static(;COMMON...,
    maximum_aerodynamic_time=.2,minimum_aerodynamic_time=.05,aerodynamic_window_time=.01,
    aerodynamic_force_tolerance=1e-4,aerodynamic_circulation_tolerance=1e-4,
    maximum_static_iterations=100,static_relaxation=.3,geometry_tolerance=1e-5,
    load_tolerance=1e-4,consecutive_equilibrium_iterations=2,
    newton_maximum_iterations=50,newton_absolute_tolerance=1e-8,
    newton_relative_tolerance=1e-8,newton_display_iterations=false)
@assert equilibrium.converged
REFINE_ONLY || serialize(joinpath(OUT,"equilibrium.jls"),equilibrium)
original_x = copy(equilibrium.structural_problem.x)
original_gamma = copy(equilibrium.aerodynamic_problem.state.system.Γ)
function run_case(equilibrium,name,force,tolerance)
    history = Float64[]
    function check_step(i,t,structure,aero)
        # Store physical speed independently of the plotted section strain.
        push!(history,maximum(norm(e.states.V) for e in structure.model.elements))
    end
    result = run_pazy_wing_uvlm(;COMMON...,initial_equilibrium=equilibrium,
        sideslip=0.,initial_airspeed_fraction=1.,airspeed_ramp_duration=0.,duration=.5,
        settling_time=.01,perturbation_amplitude=force,perturbation_duration=.04,
        save_uvlm_history=false,uvlm_save_frequency=1,
        newton_maximum_iterations=30,newton_absolute_tolerance=tolerance==1e-4 ? 1e-5 : tolerance==1e-7 ? 1e-9 : 1e-11,
        newton_relative_tolerance=tolerance==1e-4 ? 1e-5 : tolerance==1e-7 ? 1e-9 : 1e-11,
        newton_display_iterations=false,newton_always_update_jacobian=true,
        coupling_scheme=:strong,coupling_maximum_iterations=50,coupling_relaxation=.5,
        coupling_geometry_tolerance=tolerance,coupling_load_tolerance=tolerance,
        coupling_display_iterations=false,animation_frames=2,
        save_animation_history=false,save_structural_history=false,
        save_aerodynamic_load_history=false,accepted_step_callback=check_step)
    strain = 1e6 .* (result.root_axial_strain .+ .001125 .* result.root_curvature_2)
    data = hcat(result.time,strain,result.tip_out_of_plane,result.tip_twist_degrees,
        result.coupling_iterations,result.coupling_geometry_residual,
        result.coupling_load_residual,vcat(0.,history))
    writedlm(joinpath(OUT,name*".csv"),vcat(
        ["time_s" "strain_ue" "tip_m" "twist_deg" "iterations" "geometry_residual" "load_residual" "max_section_speed_mps"],data),',')
    println("AUDIT completed ",name," mean FSI=",mean(result.coupling_iterations[2:end]))
    flush(stdout)
    return (;time=result.time,strain,tip=result.tip_out_of_plane,iterations=result.coupling_iterations)
end
open(joinpath(OUT,REFINE_ONLY ? "refinement_settings.txt" : "settings.txt"),"w") do io
    println(io,COMMON)
    println(io,"Same copied equilibrium; duration .5 s; pulse .01-.05 s; coefficients 0, 1e-5, .01 N.")
    println(io,"Default FSI tolerances 1e-4, Newton 1e-5; tight FSI 1e-7, Newton 1e-9.")
    println(io,"Refinement-only mode: $REFINE_ONLY; refined FSI 1e-9, Newton 1e-11; reuse identical serialized equilibrium.")
end
summaries = []
for (tag,tol) in (REFINE_ONLY ? (("refined",1e-9),) : (("default",1e-4),("tight",1e-7)))
    quiet = run_case(equilibrium,tag*"_quiet",0.,tol)
    for (force_tag,force) in (("small",1e-5),("larger",.01))
        forced = run_case(equilibrium,tag*"_"*force_tag,force,tol)
        delta = forced.strain .- quiet.strain
        for (lo,hi) in ((.06,.16),(.2,.3),(.39,.49))
            idx = findall(t -> lo <= t < hi,forced.time)
            rms = sqrt(mean(abs2,delta[idx]))
            push!(summaries,Any[tag,force,lo,hi,rms,rms/force,mean(forced.iterations[idx])])
        end
        @assert equilibrium.structural_problem.x == original_x
        @assert equilibrium.aerodynamic_problem.state.system.Γ == original_gamma
        writedlm(joinpath(OUT,REFINE_ONLY ? "refined_response_metrics.csv" : "paired_response_metrics.csv"),vcat(
            ["tolerance" "pulse_coefficient_N" "window_start_s" "window_end_s" "difference_RMS_ue" "RMS_per_N" "mean_iterations"],
            reduce(vcat,permutedims.(summaries))),',')
    end
end
println("SMALL SIGNAL AUDIT COMPLETED: ",OUT)
