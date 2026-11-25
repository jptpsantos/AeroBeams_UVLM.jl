import AeroBeams
import WingPropellerUVLM
using DelimitedFiles, LinearAlgebra, Printf

const scratch = ARGS[1]
const package_path = realpath(pathof(WingPropellerUVLM))
startswith(package_path,realpath(scratch)) || error("The isolated UVLM source was not loaded")
const example = normpath(joinpath(@__DIR__,"..","..","examples","pazy_aerobeams_coupled"))
include(joinpath(example,"PazyWingUVLMStaticSweep.jl"))
const output = mkpath(joinpath(@__DIR__,"output"))
BLAS.set_num_threads(1)

const speed = 10.0
const density = 1.225
const alpha = deg2rad(1.0)
const alpha_step = deg2rad(0.25)
const nc,ns = 5,20
const _,span,chord,_,_ = AeroBeams.geometrical_properties_Pazy(; mesh=:refined)
const dt = chord / (speed * nc)
const wake_parameters = (; symmetric_wing=true,maximum_wake_rows=10 * nc,dt,
    maximum_time=120 * chord / speed,minimum_time=25 * chord / speed,
    window_time=4 * chord / speed,force_tolerance=1e-5,circulation_tolerance=1e-5,
    progress_frequency=1000)

function verify_matrix_paths()
    grid = zeros(3,3,4)
    for j in 0:3,i in 0:2
        grid[:,i+1,j+1] .= (chord * i / 2,span * j / 3,0.0)
    end
    model = UVLM.create_UVLMModel(surfaces=[grid],symmetric=true,
        reference=UVLM.Reference(span*chord,chord,span,zeros(3),speed,density))
    point = UVLM.create_OperatingPoint(airspeed=speed,density=density,angleOfAttack=alpha)
    problem = UVLM.create_UVLMDynamicProblem(; model,operatingPoint=point,
        timeVector=[0.0,dt,2dt],maximumWakeRows=2,trackingTimeSteps=false)
    UVLM.begin_time_step!(problem)
    UVLM.evaluate_trial!(problem)
    UVLM.commit_time_step!(problem)
    system = problem.state.system
    fs = UVLM.freestream(point)
    for recalculate in (true,false)
        case = deepcopy(system)
        UVLM.propagate_system!(case,fs,dt; additional_velocity=nothing,repeated_points=nothing,
            nwake=case.nwake,eta=0.1,calculate_influence_matrix=recalculate,
            near_field_analysis=false,derivatives=false,advance_wake=false)
        if recalculate
            system.AIC .= case.AIC
        else
            @assert isapprox(case.AIC,system.AIC; rtol=1e-13,atol=1e-12)
        end
    end
    println("Verified agreement between complete and cached influence matrices.")
end

function run_rigid()
    verify_matrix_paths()
    reference = UVLM.Reference(span * chord,chord,span,zeros(3),speed,density)
    grid = zeros(3,nc+1,ns+1)
    for j in 0:ns,i in 0:nc
        grid[:,i+1,j+1] .= (chord * i / nc,span * j / ns,0.0)
    end
    offsets = zeros(3,nc+1,ns+1)
    weights = Matrix{Float64}(I,ns+1,ns+1)
    solver = UVLM.create_UVLMSolver(coreRadius=1e-3,wakeSheddingFraction=0.1,interaction=true)
    jobs = map((alpha-alpha_step,alpha+alpha_step)) do incidence
        Threads.@spawn begin
            operatingPoint = UVLM.create_OperatingPoint(airspeed=speed,density=density,angleOfAttack=incidence)
            wake = settle_static_uvlm(grid,offsets,weights,reference,operatingPoint,solver; wake_parameters...)
            wake.converged || error("Rigid wake failed at $(rad2deg(incidence)) degrees")
            pazy_spanwise_cn(wake.aerodynamic_problem,span,0.0)
        end
    end
    sections = map(fetch,jobs)
    slope = pazy_static_derivatives(sections,alpha_step).local_incidence.cn_alpha
    profile = slope ./ (2π)
    area = (sum(profile) - (first(profile)+last(profile))/2) / ns
    @printf("Corrected rigid wing: integrated W=%0.7f, root=%0.7f, mid=%0.7f, tip=%0.7f\n",
        area,first(profile),profile[1+ns÷2],last(profile))
    @assert 0.78 < area < 0.81
    writedlm(joinpath(output,"rigid_fixed_shape.csv"),vcat(["eta" "cn_alpha_local_per_rad" "W"],
        hcat(sections[1].eta,slope,profile)),',')
end

function run_coupled()
    static_options = (; mesh=:refined,density,chordwise_panels=nc,spanwise_panels=ns,
        symmetric_wing=true,time_step_chords=1/nc,maximum_wake_rows=10*nc,
        core_radius=1e-3,wake_shedding_fraction=0.1,aerodynamic_interaction=true,
        aerodynamic_force_tolerance=1e-5,aerodynamic_circulation_tolerance=1e-5,
        maximum_static_iterations=30,static_relaxation=1.0,geometry_tolerance=1e-4,
        load_tolerance=1e-4,load_residual_force_floor=1e-6,consecutive_equilibrium_iterations=2,
        newton_maximum_iterations=50,newton_absolute_tolerance=1e-6,
        newton_relative_tolerance=1e-6,newton_display_iterations=false,progress_frequency=1000)
    result = run_pazy_wing_uvlm_static_sweep(; airspeeds=[speed],angles_of_attack=[alpha],
        sweep_angles=[0.0],slope_modes=(:fixed_shape,),alpha_step,static_options,
        local_alpha_tolerance=1e-10,save_plots=false,output_directory=joinpath(output,"coupled"))
    case = only(result.cases)
    case.central_converged && case.fixed_shape_slope_converged || error("Coupled case failed to converge")
    slope = case.cn_alpha_local_perturbations.fixed_shape.cn_alpha
    profile = slope ./ (2π)
    area = (sum(profile) - (first(profile)+last(profile))/2) / ns
    @printf("Corrected coupled wing: integrated W=%0.7f, root=%0.7f, mid=%0.7f, tip=%0.7f\n",
        area,first(profile),profile[1+ns÷2],last(profile))
end

mode = length(ARGS) >= 2 ? ARGS[2] : "all"
mode in ("all","rigid","coupled") || error("Choose all, rigid or coupled")
mode in ("all","rigid") && run_rigid()
mode in ("all","coupled") && run_coupled()
