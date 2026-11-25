using Test
include(joinpath(@__DIR__,"..","PazyFlutterSweep.jl"))
import .PazyFlutterSweep as Flutter

@testset "Fixed-speed flutter classification" begin
    t = collect(0:0.001:2.25)
    steady = hcat(fill(0.1,length(t)),fill(3.0,length(t)))
    forces = (0.025,0.05)
    oscillation = hcat(0.02 .* sin.(2π*29.4 .* t),0.7 .* cos.(2π*29.4 .* t))
    for rate in (-0.9,0.7,0.0)
        response = oscillation .* exp.(rate .* t)
        result = Flutter.classify_equilibrium(t,steady,steady+forces[1]*response,
            steady+forces[2]*response,forces,Flutter.GROWTH_DEFAULTS)
        @test result.status == (rate < 0 ? :decay : rate > 0 ? :growth : :unresolved)
        @test result.rate ≈ rate atol=0.03
        @test result.frequency ≈ 29.4 atol=0.2
    end
    saturated = oscillation .* exp.(min.(0.7 .* t,0.8))
    result = Flutter.classify_equilibrium(t,steady,steady+forces[1]*saturated,
        steady+forces[2]*saturated,forces,Flutter.GROWTH_DEFAULTS)
    @test result.status == :unresolved
    growing = oscillation .* exp.(0.7 .* t)
    result = Flutter.classify_equilibrium(t,steady,steady+forces[1]*growing,
        steady+forces[1]*growing,forces,Flutter.GROWTH_DEFAULTS)
    @test result.status == :unresolved
    @test occursin("scale",result.reason)

    # A strong low-frequency bending transient must not dominate the
    # fitted growth envelope of the shared bending/twist mode.
    multimode = copy(growing)
    multimode[:,1] .+= 0.3 .* sin.(2π*6.2 .* t) .* exp.(-2 .* t)
    result = Flutter.classify_equilibrium(t,steady,steady+forces[1]*multimode,
        steady+forces[2]*multimode,forces,Flutter.GROWTH_DEFAULTS)
    @test result.status == :growth
    @test result.rate ≈ 0.7 atol=0.05
end

@testset "Plateaus and transition brackets" begin
    t = collect(0.001:0.001:10.0)
    steady = hcat(fill(0.15,length(t)),fill(2.0,length(t)))
    primary = sin.(2π*25 .* t)
    harmonic = cos.(2π*50 .* t)
    orbit = hcat(0.001 .* primary,0.5 .* primary+0.9 .* harmonic)
    floor = (1e-6,1e-3)
    @test Flutter.classify_lco(t,steady,floor,Flutter.LCO_DEFAULTS).status == :steady
    set = Flutter.classify_lco(t,steady+orbit,floor,Flutter.LCO_DEFAULTS)
    @test set.status == :lco
    @test set.frequency ≈ 25.0 atol=0.2
    @test Flutter.classify_lco(t,steady+orbit .* exp.(-0.12 .* t),
        floor,Flutter.LCO_DEFAULTS).status == :unresolved

    row(U,status) = (; speed_mps=U,status)
    signals = [row(43.,:decay),row(43.5,:unresolved),row(44.,:growth),row(45.,:failed),row(46.,:decay)]
    @test Flutter.observed_transitions(signals,(:decay,:growth)) == [(; first=1,last=3,has_gap=true)]
    down = [row(45.,:lco),row(44.,:unresolved),row(43.,:steady)]
    bracket = only(Flutter.lco_brackets(down,(angle_deg=5.,),(name=:seed_down,),0.5))
    @test (bracket.lower_speed_mps,bracket.upper_speed_mps) == (43.,45.)
    @test bracket.transition == :extinction
    @test bracket.status == :unresolved_between && isnan(bracket.estimate_mps)
    plan = Flutter.Stepped.stepped_plan(; airspeeds=[35.,39.5,45.],sweep_order=:downsweep,
        chord=0.0989,fixed_time_step=0.0001)
    @test getproperty.(plan.plateaus,:airspeed_mps) == [45.,39.5,35.]
    @test_throws ArgumentError Flutter.Stepped.stepped_plan(; airspeeds=[35.,35.],chord=0.0989)
    @test length(only(Flutter.run_pazy_flutter_sweep(; dry_run=true).paths)) == 4
end

@testset "Adaptive bracket and verified interpolation" begin
    options = Flutter.run_pazy_flutter_sweep(; dry_run=true).options
    function evaluate(U,variant,dt,wake)
        growth = (U-43.15)*(49.9-U)/10
        return (; angle_deg=5.0,speed_mps=U,status=growth > 0 ? :growth : :decay,
            growth_per_s=growth,frequency_hz=30.0,mode_ratio_real=0.02,mode_ratio_imag=0.03,
            trim_bending_m=0.1,trim_twist_deg=2.0,trim_strain_ue=-1200.0)
    end
    mktempdir() do root
        cases = Dict{Float64,NamedTuple}(U=>evaluate(U,:fine,options.fine_dt,6.5)
            for U in (40.0,45.0,52.0))
        bounds = NamedTuple[]
        Flutter.save_equilibrium_boundaries!(bounds,cases,evaluate,options,root)
        @test getproperty.(bounds,:transition) == [:onset,:restabilization]
        @test all(b -> b.status == :resolved && b.width_mps <= options.speed_tolerance,bounds)
        @test bounds[1].estimate_mps ≈ 43.15 atol=0.01
        @test bounds[2].estimate_mps ≈ 49.9 atol=0.01
        @test isfile(joinpath(root,"boundaries.csv"))
    end
end

if "--physical" in ARGS
    @testset "Fixed-speed cases and continued plateau callback" begin
        Flutter.BLAS.set_num_threads(1)
        staticOptions = (; maximum_aerodynamic_time=0.08,minimum_aerodynamic_time=0.02,
            aerodynamic_window_time=0.006,aerodynamic_force_tolerance=5e-4,
            aerodynamic_circulation_tolerance=5e-4,maximum_static_iterations=80,
            static_relaxation=0.5,geometry_tolerance=2e-6,load_tolerance=5e-4,
            progress_frequency=10000)
        modelOptions = (mesh=:nominal,rootTorsionCorr=false)
        mktempdir() do root
            slice = (angle_deg=0.6,pilot_speed=21.0,speeds=(20.0,20.5,22.0),
                lco_seed=21.0,lco_limits=(20.0,22.0))
            sweep = Flutter.run_pazy_flutter_sweep(; slices=(slice,),stage=:equilibrium,
                model_options=modelOptions,chordwise_panels=2,spanwise_panels=4,
                coarse_dt=0.0005,fine_dt=0.00025,wake_length_chords=2.0,
                duration=0.028,pulse_duration=0.004,pulse_start=0.0,
                growth_options=(; windows=((0.004,0.016),(0.016,0.028)),block_time=0.004),
                equilibrium_options=staticOptions,plot_results=false,
                output_directory=joinpath(root,"equilibrium"))
            @test isempty(sweep.pilots)
            @test getproperty.(sweep.equilibrium_cases,:speed_mps) == collect(slice.speeds)
            @test all(r -> r.status == :unresolved && isfinite(r.trim_bending_m),sweep.equilibrium_cases)
            case = joinpath(sweep.output_directory,"aoa0p6_V20_fine")
            @test isfile(joinpath(case,"pulse_0p025.csv"))
            seen = Float64[]
            function observe(plateau,time,values)
                push!(seen,plateau.airspeed_mps)
                @test length(time) == length(values.tip_out_of_plane_m) == 12
                @test first(time) > 0
            end
            continued = Flutter.Stepped.run_pazy_stepped_bifurcation(; airspeeds=[20.0,20.5,22.0],
                sweep_order=:upsweep,plateau_time=0.006,transient_time=0.0045,
                fixed_time_step=0.0005,wake_length_chords=2.0,
                chordwise_panels=2,spanwise_panels=4,angle_of_attack=deg2rad(0.6),
                pulse_start=0.0,pulse_duration=0.004,pulse_amplitude=0.01,
                model_options=modelOptions,equilibrium_options=staticOptions,
                responses=Flutter.LCO_RESPONSES,plateau_callback=observe,
                plot_results=false,output_directory=joinpath(root,"continued"))
            @test seen == [20.0,20.5,22.0]
            @test length(unique(getproperty.(continued.boundaries,:aerodynamic_problem_id))) == 1
            Base.eval(Flutter,:(import Plots))
            rows = [(; quantity=:equilibrium,path=:increasing_speed,transition=:onset,
                lower_speed_mps=43.0,upper_speed_mps=43.5,angle_deg=5.0,
                status=:resolved,estimate_mps=43.15)]
            Base.invokelatest(Flutter.write_boundary_plot,rows,root)
            @test filesize(joinpath(root,"pazy_flutter_boundaries.png")) > 1000
        end
    end
end
