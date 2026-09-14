# Stage 2: damping convergence starting from VERIFIED aerodynamic results.
# Physical inputs come from the selection; speed can be overridden below.
import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..", "..", "..", "..")))

module ChangConvergenceAeroelastic
include(joinpath(@__DIR__,"src","Convergence.jl"))

function aeroelastic_study()
    return (;
        selection_file = normpath(joinpath(@__DIR__,"..","..","output",
            "convergence_aerodynamic","aerodynamic_selection.toml")),
        output_directory = normpath(joinpath(@__DIR__,"..","..","output","convergence_aeroelastic")),
        dry_run = false,
        make_plots = true,

        # nothing: aerodynamic speed. Set e.g. 85.0 for this aeroelastic study only.
        # RPM scales by new_speed / aerodynamic_speed (constant advance ratio).
        speed_mps = nothing,

        # nothing: selected aerodynamic value followed by two finer values.
        # Supply explicit absolute lists to choose your own refinements.
        # Lists must start at the selected value and have >=3 ordered levels.
        # Wing span remains fixed at the selected structural element count.
        sweeps = (
            wing_chord = nothing,    # e.g. [10, 12, 15] if 10 was selected
            prop_radial = nothing,
            prop_chord = nothing,
            wake_revolutions = nothing,
            core = nothing,
            azimuth_deg = nothing,
        ),
        families = [:wing_chord,:prop_radial,:prop_chord,:wake_revolutions,:core,:azimuth_deg],

        response = (
            end_time_s = 6.0,
            trim_revolutions = 10.0,
            trim_average_revolutions = 1.0,
            impulse_start_s = nothing, # automatic: after trim_revolutions
            impulse_duration_s = 0.15,
            impulse_magnitude_nm = 1200.0,
        ),
        damping = (
            fit_start_s = nothing, # after pulse + one revolution; at least 0.75 s
            block_duration_s = 1.0,
            block_overlap = 0.9,
            frequency_band_hz = (3.0,6.0), # set for the modes being investigated
            minimum_peaks = 8,
            minimum_fit_r_squared = 0.8,
            lambda_tolerance_per_s = 0.01,
            frequency_tolerance_hz = 0.1,
        ),
        structural = (
            propeller_damping_ratio = 0.0,
            stiffness_damping_ratio = 0.0,
            damping_reference_frequency_hz = nothing,
            pitch_frequency_hz = 7.97,
            yaw_frequency_hz = 7.97,
            pitch_stiffness_nm_per_rad = 19220.0,
            yaw_stiffness_nm_per_rad = 18916.0,
            twist_frequency_hz = 12.73,
            twist_stiffness_nm_per_rad = 16835.0,
            pylon_mass_per_length_kgpm = 0.0506 * 14.5939029372064 / 0.3048,
            blade_mass_kg = 1.44,
        ),
        integration = (rho_inf=1.0, state_norm_limit=1e3, propeller_angle_limit_deg=Inf),
        coupling = (maximum_iterations=10,state_tolerance=1e-5,load_tolerance=1e-2,
            equilibrium_tolerance=1e-10,coupled_equilibrium_tolerance=1e-4,relaxation=1.0),
        propeller_moment_projection = :exact_virtual_work,
        hub_load_arm_factor = 0.5,
    )
end

function operating_selection(selection,s)
    speed=get(s,:speed_mps,nothing)
    isnothing(speed) && return selection
    speed isa Real && !(speed isa Bool) && isfinite(speed) && speed>0 ||
        throw(ArgumentError("Aeroelastic speed_mps must be finite and positive, or nothing"))
    p=selection.physical
    rpm=p.rpm*(speed/p.speed_mps)
    isfinite(rpm) && rpm>0 || throw(ArgumentError("Scaled aeroelastic RPM must be finite and positive"))
    # The verification seal belongs to the original operating point only.
    return (;physical=merge(p,(;speed_mps=Float64(speed),rpm)),
        core_mode=selection.core_mode,selected=selection.selected,
        aerodynamic_reference=selection)
end

"""Run damping convergence with an optional speed override after verifying the source selection."""
function run_aeroelastic(s=aeroelastic_study())
    reference=Convergence.load_selection(s.selection_file)
    selection=operating_selection(reference,s)
    println("Aeroelastic operating point: $(selection.physical.speed_mps) m/s, $(selection.physical.rpm) RPM")
    Convergence.load_model()
    result=Base.invokelatest(Convergence.run_aeroelastic_loaded,s,selection)
    description="Aerodynamic reference: $(reference.physical.speed_mps) m/s, $(reference.physical.rpm) RPM.\n\nAeroelastic operating point: $(selection.physical.speed_mps) m/s, $(selection.physical.rpm) RPM.\n\n"
    if selection.physical.speed_mps!=reference.physical.speed_mps
        description*="The aerodynamic selection was verified at the reference speed; this study checks aeroelastic sensitivity at the overridden speed with proportional RPM.\n\n"
    end
    open(joinpath(result.directory,"report.md"),"a") do io
        print(io,"\n",description)
    end
    return merge(result,(;aerodynamic_reference=reference))
end

main() = run_aeroelastic()
end

if abspath(PROGRAM_FILE) == @__FILE__
    ChangConvergenceAeroelastic.main()
end
