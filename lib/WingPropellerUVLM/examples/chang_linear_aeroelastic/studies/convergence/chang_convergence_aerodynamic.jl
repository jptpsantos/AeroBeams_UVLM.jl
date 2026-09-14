# Stage 1: independent rigid wing/propeller convergence using CT and CQ.
# Edit parameters HERE. This script does not read chang_case.jl or CHANG_* ENV.
import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..", "..", "..", "..")))

module ChangConvergenceAerodynamic
include(joinpath(@__DIR__,"src","Convergence.jl"))

function aerodynamic_study()
    speed_mps = 65.0
    physical = (
        speed_mps = speed_mps,
        density_kgpm3 = 1.225,
        alpha_deg = 3.0,
        beta_deg = 0.0,
        wing_span_m = 7.5,
        wing_root_chord_m = 1.8,
        wing_tip_chord_m = 1.8,
        wing_symmetric = true, # Wing/wake images only. Blades never use symmetry.
        elastic_axis_fraction = 0.30,
        propeller_radius_m = 1.15,
        propeller_chord_m = 0.197,
        blades = 4,
        attachment_eta = 0.83,
        pylon_length_m = 5.6 * 0.3048,
        rpm = 1212.0 * speed_mps / 65.0, # Actual RPM at speed_mps; not a reference RPM.
        collective_offset_deg = 0.0,
        interaction = false,
        wake_shedding_fraction = 0.1,
    )

    # Core mode: :fixed (metres), :chord (fraction of full local chord),
    # or :segment (fraction of the local 3D span/radial vortex edge).
    core_mode = :fixed
    nominal = (wing_span=30, wing_chord=10, prop_radial=10, prop_chord=10,
        wake_revolutions=2.0, core=0.001, azimuth_deg=5.0)

    # Absolute panel counts per direction (propeller counts are per blade).
    # Each family varies one parameter; others retain the nominal values.
    # Use at least three ordered levels for every enabled family.
    sweeps = (
        wing_span = [20, 40, 60],
        wing_chord = [5, 10, 20],
        prop_radial = [5, 10, 20],
        prop_chord = [5, 10, 20],
        wake_revolutions = [1.0, 2.0, 3.0],
        core = [0.01, 0.003, 0.001, 0.0003],
        azimuth_deg = [5.0, 2.5, 1.25],
    )
    return (;
        physical, core_mode, nominal, sweeps,
        families = collect(keys(sweeps)),
        simulated_revolutions = 8,
        averaged_revolutions = 2,
        absolute_tolerances = (CT=1e-6, CQ=1e-6),
        relative_tolerance = 0.01,
        periodic_tolerances = (CT=1e-6, CQ=1e-6),
        dry_run = false, # true writes only the case matrix.
        make_plots = true,
        reuse_existing = true,
        # Publish a selection only after every family AND a combined
        # candidate/finer verification satisfy the CT/CQ criteria.
        confirm_selection = true,
        output_directory = normpath(joinpath(@__DIR__,"..","..","output","convergence_aerodynamic")),
    )
end

main() = Convergence.run_aerodynamic(aerodynamic_study())
end

#if abspath(PROGRAM_FILE) == @__FILE__
    ChangConvergenceAerodynamic.main()
#end
