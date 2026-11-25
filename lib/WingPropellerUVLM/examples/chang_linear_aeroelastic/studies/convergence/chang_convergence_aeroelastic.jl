# Aeroelastic convergence. Edit camelCase settings in create_study().
# Activate the local studies project if needed when loading from the IDE/REPL.
studiesProject = normpath(joinpath(@__DIR__,"..","..","..","..","studies"))
if Base.active_project() != joinpath(studiesProject,"Project.toml") &&
   (abspath(PROGRAM_FILE) == (@__FILE__) || Base.find_package("WingPropellerUVLMStudies") === nothing)
    import Pkg
    Pkg.activate(studiesProject)
end

module ChangConvergenceAeroelastic
import WingPropellerUVLMStudies
const Convergence = WingPropellerUVLMStudies.Convergence
const MovingBlockStudy = WingPropellerUVLMStudies.MovingBlockStudy
const MovingBlockDamping = WingPropellerUVLMStudies.MovingBlockDamping
operating_selection(selection,s) = WingPropellerUVLMStudies.operating_selection(selection,s)
include(joinpath(@__DIR__,"src","ExplicitAeroelasticStudy.jl"))
include(joinpath(@__DIR__,"src","StagedAeroelasticStudy.jl"))
function run_aeroelastic(s=aeroelastic_study())
    validate_run_mode(s)
    return isnothing(s.selection_file) ?
           run_explicit_aeroelastic(s) : WingPropellerUVLMStudies.run_aeroelastic(s)
end

# Select :wing for tip bending/twist, or :propeller for propeller pitch/yaw.
function create_study(; runMode=:full,dryRun=false,dampingTarget=:wing,wallHours=36.0)
    # :preview checks the coupled solver on a tiny mesh and time interval.
    # :full uses the convergence controls below and estimates damping.
    runMode in (:preview,:full,:staged) || throw(ArgumentError("runMode must be :preview, :full or :staged"))
    Convergence.Damping.response_channels(dampingTarget)
    airspeed = 65.0
    trimSpeed = 65.0
    trimRpm = 1232.63301527
    physical = (
        airspeed=airspeed,
        densityKgpm3=1.225,
        alphaDeg=3.0,
        betaDeg=0.0,
        wingSpanM=7.5,
        wingRootChordM=1.8,
        wingTipChordM=1.8,
        wingSymmetric=false,
        elasticAxisFraction=0.30,
        propellerRadiusM=1.15,
        propellerChordM=0.197,
        blades=4,
        attachmentEta=0.83,
        pylonLengthM=5.6 * 0.3048,
        rpm=trimRpm * airspeed / trimSpeed,# Constant advance ratio.
        collectiveOffsetDeg=0.0,
        interaction=false,
        wakeSheddingFraction=0.1
    )
    coreMode = :fixed # :fixed (metres), :chord or :segment (fractions).

    # Fixed values for parameters outside the active sweep.
    # Explicit mode requires absolute values for every nominal control.
    # With selectionFile set, nothing instead inherits the verified selection.
    nominal = (
        wingSpanPanels=20,
        wingChordPanels=10,
        propellerRadialPanels=10,
        propellerChordPanels=10,
        wakeRevolutions=1.0,
        core=0.01 * 0.197,
        azimuthStepDeg=10.0
    )

    # Each family varies one parameter; others retain their nominal values.
    # Add/comment entries to enable/disable families. Explicit lists need at
    # least three ordered levels and may start below or above the nominal.
    # nothing generates three refinement levels from the resolved nominal.
    # wingSpanPanels also refines the structural mesh in explicit mode.
    sweeps = (
        #wingSpanPanels = [5, 10, 15, 20, 25, 30,35, 40, 50, 60, 70, 80],
        #wingChordPanels = [5, 10, 20, 30, 40, 50, 60, 70, 80],
        #propellerRadialPanels = [1, 2, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80],
        #propellerChordPanels = [1, 2, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80],
        #wakeRevolutions = [1.0, 2.0, 4.0, 6.0],
        # Explicit coupling also needs a time-step convergence check.
        azimuthStepDeg=[30.0, 20.0, 1.25],
    )

    settings = (;
        runMode,
        # nothing: use the explicit inputs above; aerodynamic convergence is
        # unverified. A TOML path enables the strict verified-selection workflow.
        selectionFile=nothing,
        physical,coreMode,
        outputDirectory=normpath(joinpath(@__DIR__,"..","..","output",
            dampingTarget == :wing ? "convergence_aeroelastic_wing" : "convergence_aeroelastic_explicit")),
        dryRun,
        makePlots=true,

        # nothing: aerodynamic speed. Set e.g. 85.0 for this aeroelastic study only.
        # RPM scales by new_speed / aerodynamic_speed (constant advance ratio).
        airspeed=airspeed,nominal,sweeps,
        families=collect(keys(sweeps)),
        response=(
            finalTime=3.0,# Eight peaks need more than the previous 2 s record.
            baselineRevolutions=5.0,
            baselineAverageRevolutions=1.0,
            impulseStartS=0.1,# nothing: after baselineRevolutions
            impulseDurationS=0.15,
            impulseMagnitudeNm=1200.0
        ),
        damping=(
            target=dampingTarget,
            # :wing: tip_displacement_m (bending), tip_twist_deg (torsion).
            # :propeller: propeller_1_pitch_deg, propeller_1_yaw_deg.
            fitStartS=nothing,# after pulse + one revolution; at least 0.75 s
            # moving_block_2.jl: rectangular, mean-removed FFT blocks shifted
            # one sample; lambda is the slope of log amplitude against time.
            blockSizeMode=:record_fraction,# or :duration for a fixed window
            blockSize=512,# doubled/halved to satisfy the ratios below
            sizeRatioLb=0.25,
            sizeRatioUb=0.50,
            peakFromStart=1,
            peakFromEnd=0,
            blockDurationS=1.0,# used only with :duration
            blockOverlap=0.9,# used only with :duration
            frequencyBandHz=(3.0,6.0),# set to the coupled mode tracked in the selected responses
            minimumPeaks=4,
            minimumFitRSquared=0.8,
            lambdaTolerancePerS=0.01,
            frequencyToleranceHz=0.1
        ),
        structural=(
            propellerDampingRatio=0.0,
            stiffnessDampingRatio=0.0,
            dampingReferenceFrequencyHz=nothing,
            pitchFrequencyHz=7.97,
            yawFrequencyHz=7.97,
            pitchStiffnessNmPerRad=19220.0,
            yawStiffnessNmPerRad=18916.0,
            twistFrequencyHz=12.73,
            twistStiffnessNmPerRad=16835.0,
            pylonMassPerLengthKgpm=0.0506 * 14.5939029372064 / 0.3048,
            bladeMassKg=1.44
        ),
        integration=(timeIntegrator=:newmark_beta,newmarkAlpha=0.05,
            rhoInf=1.0,stateNormLimit=1e3,propellerAngleLimitDeg=Inf),
        # One UVLM/structural exchange per step. Newmark-beta advances the
        # structure; maximumIterations applies only to implicit coupling.
        coupling=(scheme=:loose_explicit,maximumIterations=10,
            stateTolerance=1e-5,loadTolerance=1e-2,equilibriumTolerance=1e-10,
            coupledEquilibriumTolerance=1e-4,relaxation=1.0,verbose=false),
        propellerMomentProjection=:exact_virtual_work
    )
    if runMode == :staged
        stagedSweeps = (
            azimuthStepDeg=[10.0,5.0,2.5,1.25],
            wakeRevolutions=[1.0,2.0,3.0],
            propellerChordPanels=[5,10,15],
            propellerRadialPanels=[10,20,25],
            wingSpanPanels=[20,30,40],
            wingChordPanels=[10,15,20]
        )
        settings = merge(settings,(;
            airspeed=85.0,
            outputDirectory=normpath(joinpath(@__DIR__,"..","..","output","convergence_aeroelastic_staged")),
            makePlots=false,wallHours,
            nominal=merge(nominal,(; azimuthStepDeg=5.0,core=0.0125 * 0.197)),
            sweeps=stagedSweeps,families=collect(keys(stagedSweeps)),
            response=merge(settings.response,(; finalTime=4.0,impulseStartS=nothing))
        ))
    elseif runMode == :preview
        previewSweeps = (propellerRadialPanels=[2,3,4],)
        settings = merge(
            settings,
            (;
                outputDirectory=settings.outputDirectory * "_preview",
                makePlots=false,
                nominal=merge(
                    nominal,
                    (;
                        wingSpanPanels=4,wingChordPanels=2,
                        propellerRadialPanels=3,propellerChordPanels=1,
                        wakeRevolutions=0.5,azimuthStepDeg=10.0)
                ),
                sweeps=previewSweeps,
                families=collect(keys(previewSweeps)),
                response=merge(
                    settings.response,
                    (;
                        finalTime=0.15,baselineRevolutions=1.0,
                        baselineAverageRevolutions=1.0,impulseStartS=nothing,
                        impulseDurationS=0.02)
                ),
                damping=merge(settings.damping,(; fitStartS=0.10))
            )
        )
    end
    return create_study(settings)
end

# Legacy settings access retains the original schema for existing scripts.
aeroelastic_study(; kwargs...) = WingPropellerUVLMStudies.legacy_settings(create_study(; kwargs...))
main(; kwargs...) = WingPropellerUVLMStudies.solve!(create_study(; kwargs...)).results

function cli_main(args=ARGS)
    runMode = nothing
    dryRun = false
    wallHours = 36.0
    dampingTarget = nothing # Without a flag, use the editable create_study() default.
    for arg in args
        if arg == "--full"
            runMode = :full
        elseif arg == "--preview"
            runMode = :preview
        elseif arg == "--staged"
            runMode = :staged
        elseif arg == "--dry-run"
            dryRun = true
        elseif startswith(arg,"--wall-hours=")
            wallHours = parse(Float64,last(split(arg,'='; limit=2)))
        elseif startswith(arg,"--damping-target=")
            dampingTarget = Symbol(last(split(arg,'='; limit=2)))
        else
            throw(ArgumentError("Unknown argument $arg; use --staged, --preview, --full, --dry-run, --wall-hours=H or --damping-target=wing|propeller"))
        end
    end
    mode = something(runMode,isnothing(dampingTarget) ? :staged : :full)
    mode == :staged && !isnothing(dampingTarget) &&
        throw(ArgumentError("Staged mode jointly fits wing-tip twist and pylon pitch. Use --full to select a moving-block damping target."))
    return isnothing(dampingTarget) ? main(; runMode=mode,dryRun,wallHours) :
           main(; runMode=mode,dryRun,dampingTarget,wallHours)
end
end

if abspath(PROGRAM_FILE) == @__FILE__
    Base.invokelatest(ChangConvergenceAeroelastic.cli_main)
end
