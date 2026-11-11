# Aerodynamic convergence. Edit camelCase settings in create_study().
# Including this file defines the entry only. Use the studies project in the REPL.
studiesProject = normpath(joinpath(@__DIR__,"..","..","..","..","studies"))
if Base.find_package("WingPropellerUVLMStudies") === nothing
    import Pkg
    Pkg.activate(studiesProject)
end

module ChangConvergenceAerodynamic
import WingPropellerUVLMStudies
const Convergence = WingPropellerUVLMStudies.Convergence

function create_study()
    airspeed = 65.0
    trimSpeed = 65.0
    trimRpm = 1232.63301527
    propellerRadius = 1.15
    muProp = trimSpeed / (trimRpm * 2 * pi / 60.0) / propellerRadius
    omega = airspeed / propellerRadius / muProp
    rpm = omega * 60.0 / (2 * pi)
    physical = (
        airspeed=airspeed,
        densityKgpm3=1.225,
        alphaDeg=3.0,
        betaDeg=0.0,
        wingSpanM=7.5,
        wingRootChordM=1.8,
        wingTipChordM=1.8,
        wingSymmetric=false,# Wing/wake images only. Blades never use symmetry.
        elasticAxisFraction=0.30,
        propellerRadiusM=propellerRadius,
        propellerChordM=0.197,
        blades=4,
        attachmentEta=0.83,
        pylonLengthM=5.6 * 0.3048,
        rpm=rpm,# Actual RPM at airspeed; derived from the constant mu_prop.
        collectiveOffsetDeg=0.0,
        interaction=false,
        wakeSheddingFraction=0.1
    )

    # Core mode: :fixed (metres), :chord (fraction of full local chord),
    # or :segment (fraction of the local 3D span/radial vortex edge).
    coreMode = :fixed
    nominal = (wingSpanPanels=20,wingChordPanels=10,propellerRadialPanels=10,propellerChordPanels=10,
        wakeRevolutions=1.0,core=0.014 * 0.197,azimuthStepDeg=5)

    # Absolute panel counts per direction (propeller counts are per blade).
    # Each family varies one parameter; others retain the nominal values.
    # Use at least three ordered levels for every enabled family.
    sweeps = (
        #wingSpanPanels = [10, 20, 30],
        #wingChordPanels = [5, 10, 15],
        propellerRadialPanels=[10,20,30,40,60],
        propellerChordPanels=[10,20,30,40,60]
        #        wakeRevolutions = [2.0, 4.0, 6.0],
        #        azimuthStepDeg = [10.0, 5.0, 2.5],
    )
    settings = (;
        physical,coreMode,nominal,sweeps,
        families=collect(keys(sweeps)),
        simulatedRevolutions=9,
        averagedRevolutions=2,
        absoluteTolerances=(CT=1e-6,CQ=1e-6),
        relativeTolerance=0.01,
        periodicTolerances=(CT=1e-6,CQ=1e-6),
        dryRun=false,# true writes only the case matrix.
        makePlots=true,
        reuseExisting=true,
        # Publish a selection only after every family AND a combined
        # candidate/finer verification satisfy the CT/CQ criteria.
        confirmSelection=true,
        outputDirectory=normpath(joinpath(@__DIR__,"..","..","output","convergence_aerodynamic"))
    )
    return WingPropellerUVLMStudies.create_AerodynamicConvergenceStudy(settings)
end

# Legacy settings access retains the original schema for existing scripts.
aerodynamic_study() = WingPropellerUVLMStudies.legacy_settings(create_study())
main() = WingPropellerUVLMStudies.solve!(create_study()).results
end


if abspath(PROGRAM_FILE) == @__FILE__
    ChangConvergenceAerodynamic.main()
end

