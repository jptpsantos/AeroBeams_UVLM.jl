# Conservative lumped-inertia validation

This validation records the study that selected the control-volume spatial-
inertia remap now used by `src/chang_model_parameters.jl`. It is kept as
supporting evidence, not as a second production implementation.

The method treats stiffness as a distributed field and remaps each complete
6-by-6 nodal spatial-inertia block over target-node control volumes. This
preserves total mass, first moments, and inertia while keeping every free-node
block positive definite.

From the repository root, regenerate the comparison tables and report with:

```powershell
julia --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/validation/conservative_lumped_modal/run_modal_comparison.jl
```

Run its focused regression tests with:

```powershell
julia --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/validation/conservative_lumped_modal/runtests.jl
```

`output/` retains only the compact reference tables and report used to justify
the production choice. The retained trim/impulse summary is historical and is
labelled with the current equivalent entry point. New smoke histories and
plots are ignored.
