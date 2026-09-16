# Aeroelastic convergence. Edit camelCase settings in create_study().
# Including this file defines the entry only. Use the studies project in the REPL.
if abspath(PROGRAM_FILE) == @__FILE__
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__,"..","..","..","..","studies")))
end

module ChangConvergenceAeroelastic
import WingPropellerUVLMStudies
const Convergence=WingPropellerUVLMStudies.Convergence
const MovingBlockStudy=WingPropellerUVLMStudies.MovingBlockStudy
const MovingBlockDamping=WingPropellerUVLMStudies.MovingBlockDamping
operating_selection(selection,s)=WingPropellerUVLMStudies.operating_selection(selection,s)
run_aeroelastic(s=aeroelastic_study())=WingPropellerUVLMStudies.run_aeroelastic(s)

function create_study()
    settings = (;
        selectionFile = normpath(joinpath(@__DIR__,"..","..","output",
            "convergence_aerodynamic","aerodynamic_selection.toml")),
        outputDirectory = normpath(joinpath(@__DIR__,"..","..","output","convergence_aeroelastic")),
        dryRun = false,
        makePlots = true,

        # nothing: aerodynamic speed. Set e.g. 85.0 for this aeroelastic study only.
        # RPM scales by new_speed / aerodynamic_speed (constant advance ratio).
        airspeed = 80,

        # nothing: selected aerodynamic value followed by two finer values.
        # Supply explicit absolute lists to choose your own refinements.
        # Lists must start at the selected value and have >=3 ordered levels.
        # Wing span remains fixed at the selected structural element count.
        sweeps = (
            wingChordPanels = nothing,    # e.g. [10, 12, 15] if 10 was selected
            propellerRadialPanels = nothing,
            propellerChordPanels = nothing,
            wakeRevolutions = nothing,
            core = nothing,
            azimuthStepDeg = nothing,
        ),
        families = [:wingChordPanels,:propellerRadialPanels,:propellerChordPanels,:wakeRevolutions,:core,:azimuthStepDeg],

        response = (
            finalTime = 6.0,
            baselineRevolutions = 10.0,
            baselineAverageRevolutions = 1.0,
            impulseStartS = nothing, # automatic: after trim_revolutions
            impulseDurationS = 0.15,
            impulseMagnitudeNm = 1200.0,
        ),
        damping = (
            fitStartS = nothing, # after pulse + one revolution; at least 0.75 s
            # moving_block_2.jl: rectangular, mean-removed FFT blocks shifted
            # one sample; lambda is the slope of log amplitude against time.
            blockSizeMode = :record_fraction, # or :duration for a fixed window
            blockSize = 512, # doubled/halved to satisfy the ratios below
            sizeRatioLb = 0.25,
            sizeRatioUb = 0.50,
            peakFromStart = 1,
            peakFromEnd = 0,
            blockDurationS = 1.0, # used only with :duration
            blockOverlap = 0.9,    # used only with :duration
            frequencyBandHz = (3.0,6.0), # acceptance check on GLOBAL dominant frequency
            minimumPeaks = 8,
            minimumFitRSquared = 0.8,
            lambdaTolerancePerS = 0.01,
            frequencyToleranceHz = 0.1,
        ),
        structural = (
            propellerDampingRatio = 0.0,
            stiffnessDampingRatio = 0.0,
            dampingReferenceFrequencyHz = nothing,
            pitchFrequencyHz = 7.97,
            yawFrequencyHz = 7.97,
            pitchStiffnessNmPerRad = 19220.0,
            yawStiffnessNmPerRad = 18916.0,
            twistFrequencyHz = 12.73,
            twistStiffnessNmPerRad = 16835.0,
            pylonMassPerLengthKgpm = 0.0506 * 14.5939029372064 / 0.3048,
            bladeMassKg = 1.44,
        ),
        integration = (timeIntegrator=:newmark_beta,newmarkAlpha=0.05,
            rhoInf=1.0,stateNormLimit=1e3,propellerAngleLimitDeg=Inf),
        coupling = (scheme=:implicit_predictor_corrector,maximumIterations=10,
            stateTolerance=1e-5,loadTolerance=1e-2,equilibriumTolerance=1e-10,
            coupledEquilibriumTolerance=1e-4,relaxation=1.0,verbose=false),
        propellerMomentProjection = :exact_virtual_work,
        hubLoadArmFactor = 0.5,
    )
    return WingPropellerUVLMStudies.create_AeroelasticConvergenceStudy(settings)
end

# Legacy settings access retains the original schema for existing scripts.
aeroelastic_study()=WingPropellerUVLMStudies.legacy_settings(create_study())
main()=WingPropellerUVLMStudies.solve!(create_study()).results
end

if abspath(PROGRAM_FILE) == @__FILE__
    ChangConvergenceAeroelastic.main()
end
