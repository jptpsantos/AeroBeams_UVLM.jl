# Execute the production Chang driver with a diagnostic-only case configuration.
# The production source files are read but not edited.

const DIAGNOSTIC_DIR = @__DIR__
const EXAMPLE_DIR = normpath(joinpath(DIAGNOSTIC_DIR, ".."))
const RUN_FILE = joinpath(EXAMPLE_DIR, "run_chang_linear_aeroelastic.jl")
const CASE_FILE = joinpath(DIAGNOSTIC_DIR, "chang_diagnostic_case.jl")
const COUPLING_FILE = joinpath(DIAGNOSTIC_DIR, "chang_uvlm_coupling_diagnostic.jl")

source = read(RUN_FILE, String)
source = replace(
    source,
    "include(joinpath(@__DIR__, \"chang_case.jl\"))" =>
        "Base.include(@__MODULE__, $(repr(CASE_FILE)))",
    "include(joinpath(@__DIR__, \"chang_uvlm_coupling.jl\"))" =>
        "Base.include(@__MODULE__, $(repr(COUPLING_FILE)))",
    "const PLOT_RESULTS = true" => "const PLOT_RESULTS = false",
    "const ANIMATE_WAKE = true" => "const ANIMATE_WAKE = false",
    "maximum_wake_rows_wing=10 * nc_wing" =>
        "maximum_wake_rows_wing=TEST_WAKE_ROWS_WING",
    "maximum_wake_rows_propeller=72" =>
        "maximum_wake_rows_propeller=TEST_WAKE_ROWS_PROPELLER",
    "FCORE = (c, Δs) -> 0.5 * Δs" =>
        "FCORE = (c, Δs) -> max(TEST_FCORE_SEGMENT_FACTOR * Δs, TEST_FCORE_CHORD_FACTOR * c)",
)

# Accept either the reference half-arm or a temporary arm value in the
# production driver while keeping the diagnostic choice explicit.
source = replace(
    source,
    r"hub_center_load_A\s*=\s*SVector\(-[0-9.]+\s*\*\s*L_pylon,\s*0.0,\s*0.0\)" =>
        "hub_center_load_A = SVector(-TEST_HUB_LOAD_ARM_FACTOR * L_pylon, 0.0, 0.0)",
)

Base.include_string(@__MODULE__, source, RUN_FILE)
