# Full convergence: mesh, wake length, core radius, and time step.
# Edit this study, then run run_chang_core_convergence.jl.
# Physical conditions come from chang_case.jl; these settings only define the study.
function core_convergence_defaults(config)
    nominal = something(config.aerodynamic.core_radius_m, 1e-3)
    return (;
        mode = :aerodynamic,       # :aerodynamic or :aeroelastic
        dry_run = false,          # false runs simulations; true saves only the case matrix.
        make_plots = true,
        reuse_existing = true,    # Verified rigid results only; damping cases run fresh.
        output_directory = nothing, # nothing creates a timestamped output directory.

        # Radii in metres, largest to smallest. Include the active case value.
        radii_m = sort(unique([nominal, 0.01, 0.003, 0.001, 0.0003]); rev = true),
        nominal_radius_m = nominal,
        # Absolute panel counts, independent of the mesh in chang_case.jl.
        # Baseline used for the wake/time-step studies and unchanged directions.
        nominal_mesh = (wing_span=20, wing_chord=5, prop_radial=5, prop_chord=5),

        # Each separate mesh family changes only its corresponding direction.
        # Lists may have different lengths; use at least three increasing counts.
        panel_levels = (
            wing_span = [20, 40, 60],
            wing_chord = [5, 10, 20],
            prop_radial = [5, 10, 20],
            prop_chord = [5, 10, 20],
        ),

        # Explicit combined meshes: every radius is tested on every row.
        # These rows are independent of panel_levels above.
        core_meshes = [
            (wing_span=20, wing_chord=5,  prop_radial=5,  prop_chord=5),
            (wing_span=40, wing_chord=10, prop_radial=10, prop_chord=10),
            (wing_span=60, wing_chord=20, prop_radial=20, prop_chord=20),
        ],
        # All families enabled for the full aerodynamic convergence study.
        families = [:core_mesh, :wing_span, :wing_chord, :prop_radial,
                    :prop_chord, :time_step, :wake_length],
        # core_mesh repeats the complete radius sweep on each mesh.
        # The same results also quantify mesh convergence at each fixed radius.
        # Aeroelastic mode fixes wing span/structural resolution at
        # nominal_mesh.wing_span, overriding core_meshes.wing_span and omitting
        # the wing_span family. Other mesh directions remain independently editable.
        azimuth_levels_deg = [5.0, 2.5],
        wake_levels_revolutions = [1.0, 2.0, 3.0],
        nominal_wake_revolutions = 2.0, # Both wakes retain this age; rows scale with dt.

        simulated_revolutions = 6, # Rigid runs: fill wakes before averaging.
        averaged_revolutions = 2,
        coefficient_absolute_tolerances = (1e-4, 1e-6, 1e-6), # CL, CT, CQ
        coefficient_relative_tolerance = 0.01,
        periodic_tolerances = (1e-4, 1e-6, 1e-6),

        # Aeroelastic runs use the case's duration and excitation.
        fit_start_s = nothing,    # Automatic: after the pulse and at least one revolution.
        block_duration_s = 1.0,   # Physical time, held fixed as dt changes.
        block_overlap = 0.9,
        frequency_band_hz = (3.0, 6.0), # Track the same modal band in every run.
        minimum_fit_r_squared = 0.8,
        minimum_peaks = 8,
        lambda_tolerance_per_s = 0.01,
        frequency_tolerance_hz = 0.1,
    )
end
