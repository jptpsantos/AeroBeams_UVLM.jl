# Execute the production Chang driver with a diagnostic-only case configuration.
# The production source files are read but not edited.

const DIAGNOSTIC_DIR = @__DIR__
const EXAMPLE_DIR = normpath(joinpath(DIAGNOSTIC_DIR, ".."))
const RUN_FILE = joinpath(EXAMPLE_DIR, "run_chang_linear_aeroelastic.jl")
const CASE_FILE = joinpath(DIAGNOSTIC_DIR, "chang_diagnostic_case.jl")

source = read(RUN_FILE, String)
source = replace(
    source,
    "include(joinpath(@__DIR__, \"chang_case.jl\"))" =>
        "Base.include(@__MODULE__, $(repr(CASE_FILE)))",
    "const PLOT_RESULTS = true" => "const PLOT_RESULTS = false",
    "const ANIMATE_WAKE = true" => "const ANIMATE_WAKE = false",
    "maximum_wake_rows_wing=10 * nc_wing" =>
        "maximum_wake_rows_wing=TEST_WAKE_ROWS_WING",
    "maximum_wake_rows_propeller=72" =>
        "maximum_wake_rows_propeller=TEST_WAKE_ROWS_PROPELLER",
)

Base.include_string(@__MODULE__, source, RUN_FILE)
