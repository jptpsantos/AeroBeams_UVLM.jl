# Run directly from the IDE or the examples Julia project.
include(joinpath(@__DIR__,"PazyWingUVLMStaticSweep.jl"))

# The wake velocities are computed in parallel. The number of threads is fixed. When in VSCode, set "julia.NumThreads" (e.g., to 8 or auto) and restart the REPL; from a terminal, use `julia -t 8`.
Threads.nthreads() > 1 || @warn "Running with a single thread: the UVLM wake calculation will be serial"

# EDIT THESE FLIGHT CONDITIONS. Angles entered here are in DEGREES.
velocity_range = vcat(1e-2,10.0:10.0:120.0)  # m/s; every value must be > 0
# The first case approximates 0 deg. Use its actual incidence in the CSV.
angles_of_attack_deg = [1e-4, 1.0, 3.0, 5.0, 7.0, 8.5, 10.0]
sweep_angles_deg = [0.0, 10.0, 20.0] # native Pazy data: [0.0,10.0,20.0]
# One static equilibrium per condition, with prescribed alpha held constant.
# cn_alpha_local = cn / alpha_local, where local incidence follows the
# deformed chord/normal and uses the freestream (excluding induced downwash).
# This is a secant approximation [1/rad], not a perturbation derivative.
# Outputs use spanwise_panels + 1 nodal sections, including root and tip;
# endpoint cn values use half-width tributary areas.
slope_modes = (:fixed_shape,)  # :local_ratio, :fixed_shape, :equilibrium
# Add :fixed_shape and/or :equilibrium to compare with central differences.
# The perturbation CSV also divides delta cn by delta alpha_local at each station.
alpha_step_deg = 0.1
local_alpha_tolerance = 1e-10          # rad; supports the small-incidence case

# Native Pazy model. Same physical upstream-facing leading edge, streamwise
# chord and downstream wake as run_pazy_wing_uvlm.jl. Positive sweep is aft.
mesh = :nominal
useAePW3Matrices = false
withWingtipRod = true
rootTorsionCorr = true
beamWingtipRod = false
nElementsWingtipRod = 12
sweptTipMassConfig = :baseline
sweptTipMassInertiaCorrection = true
density = 1.225                         # kg/m^3

# UVLM mesh and free wake. Refine these for quantitative convergence studies.
chordwise_panels = 5
spanwise_panels = 20
symmetric_wing = true                   # false: isolated cantilever
match_wake_panel_chord = true           # nominal wake-row length equals one wing chordwise panel
time_step_chords = 0.25
if match_wake_panel_chord
    time_step_chords = 1 / chordwise_panels
end
wake_length_chords = 10.0
maximum_wake_rows = ceil(Int,wake_length_chords / time_step_chords)
core_radius = 1e-6                      # m
wake_shedding_fraction = 0.1
aerodynamic_interaction = true

# Wake convergence times in U*t/chord, scaled separately for every speed.
# Each static trial starts a fresh wake on a fixed shape. A minimum of two
# retained wake ages is required by the static solver, in addition to these.
maximum_aerodynamic_chords = 120.0
minimum_aerodynamic_chords = 25.0
aerodynamic_window_chords = 4.0
aerodynamic_force_tolerance = 1e-5
aerodynamic_circulation_tolerance = 1e-5

# Nonlinear static equilibrium (no structural inertia or dynamic settling).
maximum_static_iterations = 100
# Load update of the static iteration. :anderson (Anderson mixing, mixing parameter static_relaxation = 1) converges in a few iterations with no tuning; :fixed uses loads <- (1-w) loads + w G(loads), w = static_relaxation (0.3 is safe but slow, 1 can diverge).
relaxation_mode = :anderson
anderson_depth = 4
static_relaxation = 1.0
# Maximum UVLM grid change per load update, in chords (the step is halved above it).
geometry_cap_chords = 2.0
# Keep the wake between static iterations (and between perturbed solves) instead of a fresh one.
wake_reuse = true
# Start each case from the previous one (same sweep angle) and perturbed solves from the central one.
warm_start = true
geometry_tolerance = 1e-4
load_tolerance = 1e-4
load_residual_force_floor = 1e-6        # N; resolve the tiny loads in the zero-incidence proxy
consecutive_equilibrium_iterations = 2
newton_maximum_iterations = 50
newton_absolute_tolerance = 1e-6
newton_relative_tolerance = 1e-6
newton_display_iterations = false
progress_frequency = 1000

# Base folder; each prescribed incidence gets its own aoa...deg subfolder.
# Examples: aoa3deg, aoa2p5deg, aoam2p5deg.
output_directory = joinpath(@__DIR__,"output","static_spanwise_cn_cm")
save_plots = true

static_options = (; mesh,useAePW3Matrices,withWingtipRod,rootTorsionCorr,
    beamWingtipRod,nElementsWingtipRod,sweptTipMassConfig,sweptTipMassInertiaCorrection,
    density,chordwise_panels,spanwise_panels,symmetric_wing,time_step_chords,
    maximum_wake_rows,core_radius,wake_shedding_fraction,aerodynamic_interaction,
    aerodynamic_force_tolerance,aerodynamic_circulation_tolerance,
    maximum_static_iterations,static_relaxation,geometry_tolerance,load_tolerance,load_residual_force_floor,
    consecutive_equilibrium_iterations,newton_maximum_iterations,
    newton_absolute_tolerance,newton_relative_tolerance,newton_display_iterations,progress_frequency,
    relaxation_mode,anderson_depth,geometry_cap_chords,wake_reuse)

static_sweep = run_pazy_wing_uvlm_static_sweep(;
    airspeeds=velocity_range,angles_of_attack=deg2rad.(angles_of_attack_deg),
    sweep_angles=deg2rad.(sweep_angles_deg),
    slope_modes,alpha_step=deg2rad(alpha_step_deg),local_alpha_tolerance,static_options,
    maximum_aerodynamic_chords,minimum_aerodynamic_chords,
    aerodynamic_window_chords,output_directory,save_plots,warm_start)

println("Spanwise cn, cn_alpha and cm outputs, grouped by angle of attack under ",output_directory)
println("Central equilibria converged: ",count(case -> case.central_converged,static_sweep.cases),
    "/",length(static_sweep.cases))
# Example REPL access:
# static_sweep.cases[1].section.eta
# static_sweep.cases[1].section.alpha_local   # radians, varies along the deformed wing
# static_sweep.cases[1].cn_alpha_local       # cn / alpha_local, 1/rad
# static_sweep.cases[1].section.cm_quarter_chord
# static_sweep.cases[1].section.cm_spar
# Ratios are NaN for unconverged cases or near-zero local incidence.
# Inspect cn_alpha_local_valid in the CSV and the convergence histories.
