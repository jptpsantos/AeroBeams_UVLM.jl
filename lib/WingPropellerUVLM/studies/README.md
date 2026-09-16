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

The aerodynamic and aeroelastic entry scripts remain in their original locations. Edit their camelCase settings in `create_study()`. Including an entry defines functions; executing it directly runs the study. Legacy `aerodynamic_study()` and `aeroelastic_study()` return the old NamedTuple schema.

Set `dryRun=false` to execute a refinement. An aeroelastic study requires a verified aerodynamic selection under the current source. Its optional `airspeed` override retains the original reference and scales RPM proportionally. The shared runner uses the reference moving-block method.

The standalone aerodynamic package exposes `create_UVLMModel`, `create_OperatingPoint`, `create_UVLMSolver`, `create_UVLMDynamicProblem`, and `solve!`. For a single Chang case, edit `examples/chang_linear_aeroelastic/chang_case.jl` and run `run_chang_linear_aeroelastic.jl`. The studies API also accepts those defaults through `canonical_settings` and `create_ChangModel`. Native AeroBeams structural coupling remains separate work.

The complete self-contained LaTeX guide is [UVLM_ANALYSIS_GUIDE.tex](../../../docs/migration/UVLM_ANALYSIS_GUIDE.tex). It explains every active code area, units, settings, execution and result interpretation, and embeds the actual entry code.

Run all library and study tests in the declared environment:

```sh
julia --project=lib/WingPropellerUVLM/studies lib/WingPropellerUVLM/studies/test/runtests.jl
```
