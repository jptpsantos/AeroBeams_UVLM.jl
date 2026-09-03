# Run the production Chang aeroelastic driver with the conservative lumped
# inertia remap injected in memory. Production source files are read, not
# edited. Mesh, flow, wake, core, and duration are controlled by the same
# CHANG_TEST_* variables used by the existing diagnostic harness.

const EXPERIMENT_DIR = @__DIR__
const EXAMPLE_DIR = normpath(joinpath(EXPERIMENT_DIR, "..", ".."))
const DIAGNOSTIC_DIR = joinpath(EXAMPLE_DIR, "diagnostics")
const RUN_FILE = joinpath(EXAMPLE_DIR, "run_chang_linear_aeroelastic.jl")
const CASE_FILE = joinpath(DIAGNOSTIC_DIR, "chang_diagnostic_case.jl")
const COUPLING_FILE = joinpath(DIAGNOSTIC_DIR, "chang_uvlm_coupling_diagnostic.jl")
const REMAP_FILE = joinpath(EXPERIMENT_DIR, "apply_conservative_remap.jl")

source = read(RUN_FILE, String)
parameter_include = "include(joinpath(@__DIR__, \"chang_model_parameters.jl\"))"
source = replace(
    source,
    "include(joinpath(@__DIR__, \"chang_case.jl\"))" =>
        "Base.include(@__MODULE__, $(repr(CASE_FILE)))",
    "include(joinpath(@__DIR__, \"chang_uvlm_coupling.jl\"))" =>
        "Base.include(@__MODULE__, $(repr(COUPLING_FILE)))",
    parameter_include =>
        parameter_include * "\nBase.include(@__MODULE__, $(repr(REMAP_FILE)))",
    "const PLOT_RESULTS = true" => "const PLOT_RESULTS = false",
    "const ANIMATE_WAKE = true" => "const ANIMATE_WAKE = false",
    "maximum_wake_rows_wing=10 * nc_wing" =>
        "maximum_wake_rows_wing=TEST_WAKE_ROWS_WING",
    "maximum_wake_rows_propeller=72" =>
        "maximum_wake_rows_propeller=TEST_WAKE_ROWS_PROPELLER",
    "FCORE = (c, Δs) -> 0.5 * Δs" =>
        "FCORE = (c, Δs) -> max(TEST_FCORE_SEGMENT_FACTOR * Δs, TEST_FCORE_CHORD_FACTOR * c)",
)

source = replace(
    source,
    r"hub_center_load_A\s*=\s*SVector\(-[0-9.]+\s*\*\s*L_pylon,\s*0.0,\s*0.0\)" =>
        "hub_center_load_A = SVector(-TEST_HUB_LOAD_ARM_FACTOR * L_pylon, 0.0, 0.0)",
)

occursin(repr(REMAP_FILE), source) || error("Failed to inject conservative remap")
Base.include_string(@__MODULE__, source, RUN_FILE)
