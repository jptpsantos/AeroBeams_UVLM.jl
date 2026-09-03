# Conservative lumped-inertia modal experiment

This folder is intentionally independent of the active Chang aeroelastic
implementation. It does not modify or include `chang_model_parameters.jl`.

The model applies two distinct spanwise treatments:

- stiffness is a distributed field sampled at each beam-element center;
- mass, center of gravity, and inertia are remapped as complete 6-by-6 spatial
  inertia blocks and lumped at every free target node.

The 16 spreadsheet records and the zero root record are the nodal reference
values. Each complete nodal spatial-inertia block is divided by its source
tributary length. The resulting density is linearly interpolated between
reference stations and integrated exactly over every new beam element. The
integrated property is lumped at the element's outer node, so a uniform mesh
may use any positive number of elements. Because all interpolation weights are
nonnegative, valid source spatial inertias cannot become indefinite after the
remap. The fixed root remains massless because its six DOFs are constrained;
every free node has a positive lumped spatial inertia.

Run the analysis from the repository root:

```powershell
julia --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/experiments/conservative_lumped_modal/run_modal_comparison.jl
```

Run the regression tests with:

```powershell
julia --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/experiments/conservative_lumped_modal/runtests.jl
```

The comparison script writes its CSV tables and Markdown report under
`output/`. It also reproduces the present 20-element failure using the current
signed-component rescaling algorithm and contrasts it with the former
unscaled interpolation.

To run the actual production aeroelastic driver with this remap injected in
memory, use the separate diagnostic launcher. For example, a short 20-element
smoke test in PowerShell is:

```powershell
$env:CHANG_TEST_WING_SPAN_PANELS = "20"
$env:CHANG_TEST_WING_CHORD_PANELS = "3"
$env:CHANG_TEST_END_TIME_S = "0.01"
$env:CHANG_END_TIME_S = "0.01"
$env:CHANG_OUTPUT_DIR = "lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/experiments/conservative_lumped_modal/output/aeroelastic_20"
$env:CHANG_OUTPUT_LABEL = "conservative_lumped_20"
julia --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/experiments/conservative_lumped_modal/run_conservative_diagnostic_case.jl
```

The launcher reads the current production driver, substitutes the diagnostic
case/coupling files, disables plots, and injects `apply_conservative_remap.jl`
immediately after the normal parameter file. It does not edit the production
files.
