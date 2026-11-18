# UVLM studies

This package owns shared convergence, moving-block damping and the retained Chang benchmark interface. Start from the repository root:

```sh
julia --project=lib/WingPropellerUVLM/studies lib/WingPropellerUVLM/studies/setup.jl
julia --project=lib/WingPropellerUVLM/studies
```

```julia
import WingPropellerUVLMStudies as Studies
include("lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/studies/convergence/chang_convergence_aerodynamic.jl")
base = ChangConvergenceAerodynamic.create_study()
study = Studies.create_AerodynamicConvergenceStudy(base.settings;
    dryRun=true, makePlots=false, outputDirectory="output/check_matrix")
Studies.solve!(study)
```

The aerodynamic and aeroelastic entry scripts remain in their original locations. Edit their camelCase settings in `create_study()`. The aerodynamic preset calls `main()` at the bottom; comment that call when loading its definitions for manual REPL use. The aeroelastic entry runs only when executed directly, so including it loads definitions without starting a study. Legacy `aerodynamic_study()` and `aeroelastic_study()` return the old NamedTuple schema.

Set `dryRun=false` to execute a refinement. The Chang aeroelastic entry defaults to explicit physical and nominal inputs with `selectionFile=nothing`, and labels aerodynamic convergence as unverified. Its entry wrapper also supports a span sweep that rebuilds both aerodynamic and structural meshes. Set `selectionFile` to a TOML path to use the shared runner's strict verified-selection workflow under the current source. Its optional `airspeed` override retains the original reference and scales RPM proportionally. Both entry workflows use the reference moving-block method.

The aeroelastic entry's `dampingTarget=:wing` selects wing-tip bending displacement and twist; `:propeller` selects propeller pitch and yaw. The shared settings use `damping.target`, with `:propeller` retained when old settings omit the field. CSV metric names, diagnostics and fit metadata identify the selected responses and units. Use `--damping-target=wing` or `--damping-target=propeller` on the entry command line. See the [convergence settings and signal definitions](../examples/chang_linear_aeroelastic/studies/convergence/README.md#full-explicit-analysis-and-runtime-planning).

The standalone aerodynamic package exposes `create_UVLMModel`, `create_OperatingPoint`, `create_UVLMSolver`, `create_UVLMDynamicProblem`, and `solve!`. For a single Chang case, edit `examples/chang_linear_aeroelastic/chang_case.jl` and run `run_chang_linear_aeroelastic.jl`. The studies API also accepts those defaults through `canonical_settings` and `create_ChangModel`. Native AeroBeams structural coupling remains separate work.

The complete self-contained LaTeX guide is [UVLM_ANALYSIS_GUIDE.tex](../../../docs/migration/UVLM_ANALYSIS_GUIDE.tex). It explains every active code area, units, settings, execution and result interpretation, and embeds the actual entry code.

Run all library and study tests in the declared environment:

```sh
julia --project=lib/WingPropellerUVLM/studies lib/WingPropellerUVLM/studies/test/runtests.jl
```
