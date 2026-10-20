# Efficiency implementation — 19 October 2026

This records the changes requested after the [efficiency review](README.md). The implementation removes repeated allocation and computation in the existing solvers, simplifies duplicated orchestration, and makes detailed output optional. It does not change the physical models or the user's numerical settings.

## Scope and reference

The review traced AeroBeams structural assembly, Newton solution, distributed loads, strip aerodynamics and AD; the UVLM geometry, circulation, near-field force and wake paths; Pazy coupling and LCO drivers; and the general wing/propeller structural adapter and time integrators. The original review lists additional architectural opportunities that need separate work.

The implementation baseline is commit `0b02c08730d8072ea46383766863bc888e761017`. The short coupled reference in [coupled_baseline.toml](coupled_baseline.toml) was captured with that source before production edits. The earlier isolated measurements in [measurements.toml](measurements.toml) and [cpu_profile.txt](cpu_profile.txt) are preserved.

The current LCO driver still uses 45–60 m/s in 5 m/s increments, 3 degrees angle of attack, zero perturbation amplitude, 10 seconds per speed, the same mesh, wake length, timestep, relaxation and convergence tolerances. The full LCO sweep was not run for this implementation. An existing `output/LCO_sweep/pazy_LCO_V45.0.csv` was left untouched.

## AeroBeams changes

| Source | Modification and reason |
| --- | --- |
| [`src/Problem.jl`](../../../src/Problem.jl) | Add one shared inventory of the eight distributed load tables. Precompute only configured categories, keeping inactive tables as `(2,3,0)` arrays. Recompute element flags from the resulting values, including resetting flags when a category becomes zero. Retain adaptive quadrature for active categories. Save independent copies in `xOverTime`. |
| [`src/Element.jl`](../../../src/Element.jl) | Declare these tables as `Array{Float64,3}` so their dimensionality is concrete. |
| [`src/SystemSolver.jl`](../../../src/SystemSolver.jl) | Add `checkpoint_model`: deep-copy the mutable model graph while sharing precomputed load tables and `R_A_ofTime`, which Newton reads without mutating. Copy state/residual vectors normally and assign scalar checkpoints directly. Add a sparse Newton factorization cache local to one nonlinear solve. Use `ldiv!` into the increment and negate that vector. Reuse numeric factors only while the Jacobian is unchanged and symbolic factors only while the sparse pattern matches. Reset factors for each load step/retry. Preserve spring Jacobians when skipping a Jacobian update. Use views for state updates. Save independent copies in `xOverσ`. |
| [`src/Core.jl`](../../../src/Core.jl) | Gather states through views, fuse elementwise rate/scaling expressions, and clear owned gravity buffers. Pass the load factor to interpolation and scale only the six returned values, preserving the existing time/index convention and owned outputs. Remove redundant `copy` operations before `.=` and materialized transpose matrices before `mul3`. |
| [`src/Utilities.jl`](../../../src/Utilities.jl) | Store the six fixed Gauss–Legendre rules once instead of constructing a dictionary and every rule per call. Build `mul3` in one owned output matrix using three `mul!` calls. |
| [`src/AeroBeams.jl`](../../../src/AeroBeams.jl) | Remove eager `Revise` loading from numerical package import. Its dependency remains available to development environments. |

### Checkpoint and factorization contracts

`checkpoint_model` deliberately retains ordinary `deepcopy` behavior for mutable element, BC, beam, constraint and load-cache state. It preserves internal graph aliases, including `element === beam.elements[i]` and BC references to the copied beam. It avoids a hand-maintained inventory of every mutable Newton field. Model replacement on a retry follows the previous behavior.

The shared load arrays must remain read-only inside a nonlinear solve. `precompute_distributed_loads!` replaces them outside the solve. `R_A_ofTime` is output history and is also read-only inside Newton. Do not append to or mutate either through a checkpoint. A custom extension that mutates them during Newton must revise this contract. The implementation uses Julia's internal `Base.deepcopy_internal` to seed the identity map; this requires checking when changing supported Julia versions.

The sparse cache compares both `colptr` and `rowval` before symbolic reuse. A changed pattern creates a new factorization. This also avoids a Julia 1.11.5 stdlib error encountered with `lu!(F,A; reuse_symbolic=false)`. Trim and hinge-constrained solution branches retain their existing specialized solves. This change does not introduce a preassembled structural sparsity pattern.

Saved `xOverTime` and `xOverσ` entries now own their vectors. This is a correctness fix: the previous entries could alias the live state and change after saving. It adds the necessary state-vector copy when a history is requested.

## UVLM backend and transaction changes

Paths in this section are relative to `lib/WingPropellerUVLM/`.

| Source | Modification and reason |
| --- | --- |
| `src/backend/system.jl` | Allocate `span_force_velocities` alongside existing segment force buffers. |
| `src/backend/nearfield.jl` | Accept optional chord, span, unsteady and segment-velocity buffers. Reuse them in the System path, validate their dimensions, and clear entries that can otherwise retain stale forces. Evaluate circulation jumps from the original circulation with the equivalent sign instead of allocating its negation. Calls without buffers retain allocating behavior. |
| `src/backend/analyses.jl` | Pass reusable buffers to the default near-field routine in steady/dynamic analyses. Preserve legacy callback keywords and custom callback behavior. Remove the trailing-edge update immediately following a full AIC rebuild, since the full assembly already supplies those columns. Keep the trailing-edge-only update when reusing a fixed-geometry AIC. Detect actual surface motion instead of setting it unconditionally true; moving surfaces still force a full rebuild. |
| `src/backend/geometry.jl` | Calculate chord lengths from fixed three-component differences instead of allocating slices in panel construction/update. |
| `src/UVLMState.jl` | Add `copy_system!` to reset all history and scratch arrays into existing storage for a fixed mesh/wake capacity. Mutable nonlinear section records retain independent copies. The existing snapshot/restore API remains available. |
| `src/UVLMProblem.jl` | Allocate a separate commit System once. Reset trial storage from accepted state in place; update panels using stored geometry/ratios. Avoid the shared trailing-edge search in trials that defer wake advancement. Convection runs in the private commit buffer, and accepted/commit buffers swap only on success. Rollback resets existing trial storage. Advance wake row counts without temporary arrays. |
| `src/backend/analyses.jl`, `circulation.jl`, `induced.jl`, `wake.jl` | Remove obsolete commented function bodies. Remove ineffective duplicate-velocity assignments in wake evaluation: their inner-loop `continue` never skipped the actual velocity calculation and the copied values were overwritten. Actual shared-vertex computation is retained. This cleanup makes no claim of implementing cross-surface velocity reuse. |

### UVLM ownership and lifetime

Accepted, trial and commit storage remain isolated. Every trial starts with the same accepted circulation, geometry and wake. Repeated trials do not convect the wake. A failed geometry trial or failed wake commit leaves accepted state unchanged; a failed commit remains retryable. Successful convection occurs once per accepted physical step.

`copy_system!` assumes fixed surface count, mesh shape and allocated wake capacity. A topology/capacity change requires creating a new problem. The facade already checks trial topology before applying geometry.

Systems and scratch arrays belong to the problem. The previous accepted System becomes a reusable commit buffer, so an external reference to it is not a permanent historical snapshot. Use the owned result histories, `snapshot_uvlm`, or an explicit copy to retain data. Trial loads returned by `aerodynamic_loads` retain owned output arrays as before.

The facade now retains three Systems (accepted, trial and commit). The extra persistent commit storage trades a modest fixed memory increase for avoiding repeated full-System allocation and preserving failed-commit isolation. It does not make every backend or returned-load operation allocation-free.

The four-argument `UVLMWorkspace` constructor remains supported and creates the extra private buffer. Direct positional construction of the backend `System` must account for the new scratch field; normal `System(...)` grid/surface constructors do this automatically. No serialized backend System compatibility is claimed.

The full AIC is still assembled for every facade coupling trial. No factorization is cached across changing aerodynamic geometry. Fixed-geometry backend callers requesting matrix reuse must perform an initial full assembly and rebuild it when interaction/symmetry/core/trailing-edge policies change.

## Pazy coupling, history and drivers

| Source in `examples/pazy_aerobeams_coupled/` | Modification and reason |
| --- | --- |
| `PazyWingUVLMCoupling.jl` | Precompute `PazySpanMap`, containing each station's nonzero beam-node support and weights. Add `wing_geometry!` and `beam_loads!`, retain allocating wrappers, and reuse coupling geometry/load buffers. Use fixed three-component vectors and direct coordinate writes. Preserve the clamped root and each node's own vortex lever arm. Compute load/geometry residuals by scalar reductions. Keep the final consistency trial. Use a load column view for BC updates. |
| `PazyWingUVLMCoupling.jl` | Add `save_structural_history=true` and `save_aerodynamic_load_history=true`. Skip structural snapshots and orientation-history appends when disabled; disabled nodal load history has shape `(6,n_nodes,0)`. Scalar tip responses, time, airspeed and coupling diagnostics are still saved. Include the flags in the result. |
| `PazyWingUVLMVisualization.jl` | Report missing required histories before rendering. Permit wake-only frames without structural snapshots; synchronize structural/wake frames when structural history exists. Require structural and nodal load histories for force glyphs. Avoid scanning force history when no vectors are requested. |
| `run_pazy_wing_uvlm.jl` | Request structural and nodal load history when generating animations. Response plots still use scalar histories. |
| `run_pazy_LCO.jl` | Wrap orchestration in `run_pazy_LCO()` so sweep variables are local and each completed simulation can be released after extracting its response. Store only response summaries/histories in the returned `lco_results`. Explicitly disable detailed structural and nodal load history alongside the existing disabled UVLM/animation histories. Keep direct execution from the IDE. |

The span map stores arbitrary nonzero support, rather than assuming exactly two entries. The current linear interpolation uses at most two. Geometry interpolation and reverse load transfer use the same weights and node-specific offsets, preserving the existing force/moment and virtual-work mapping.

Detailed history defaults remain enabled for direct `run_pazy_wing_uvlm` callers. An existing caller that needs structural plots can continue using those defaults. For a response-only analysis, pass:

```julia
save_structural_history = false,
save_aerodynamic_load_history = false,
save_animation_history = false,
save_uvlm_history = false,
```

The LCO driver sets these flags itself. `lco_results` exposes `speeds`, `response_amplitudes`, `speed_histories`, `aoa_histories` and `output_directory`. Sweep internals formerly left as script globals are now local to the function.

## General wing/propeller model and linear integrators

| Source under `lib/WingPropellerUVLM/` | Modification and reason |
| --- | --- |
| `src/aeroelastic/WingPropellerStructure.jl` | Scatter local attachment mass and gyroscopic blocks directly into node DOFs instead of multiplying through dense selector matrices. Keep returned attachment operators for compatibility. Reconstruct aerodynamic kinematics with views and one owned output per component. Allow already reconstructed kinematics to be supplied to hub-position evaluation. |
| `examples/wing_propeller_aeroelastic/src/general_uvlm_coupling.jl` | Reuse the trial kinematics for propeller hubs instead of reconstructing them a second time. |
| `src/aeroelastic/GeneralizedAlpha.jl` | Delegate the older partitioned generalized-alpha entry point to the common partitioned loop with the same method, arguments and acceptance semantics. Retain integrator-specific formulas. |
| `src/aeroelastic/StructuralTimeIntegration.jl` | Accept an optional effective system in loose/partitioned steps. Build it only if the caller has not supplied one. |
| `examples/wing_propeller_aeroelastic/src/WingPropellerAeroelastic.jl` | Retain the effective structural system across physical steps for the fixed run-level `M,C,K` and integration parameters. Rebuild when `dt` changes. Use exact timestep equality, avoiding approximate reuse for genuinely different steps. |

An externally supplied effective system is valid only for the matrices, timestep and method parameters used to construct it. Callers changing any of these must rebuild it. The common driver constructs its structural model once per run; callbacks should not mutate that model's matrices.

## Verification and reproducibility

The following checks are included in the normal test entry points:

- [`test/efficiencyTests.jl`](../../../test/efficiencyTests.jl), included by the root `test/runtests.jl`: inactive history storage, analytical active distributed loads, scaled interpolation at step boundaries/interior times, output ownership, checkpoint graph/state isolation, sparse cache invalidation, saved-state ownership, all fixed quadrature rules, `mul3`, an analytical spring with modified Newton, and a selected Pazy strip-AD/tip-correction reference case.
- [`lib/WingPropellerUVLM/test/efficiency_regressions.jl`](../../../lib/WingPropellerUVLM/test/efficiency_regressions.jl), included by its `test/runtests.jl`: full AIC versus trailing-edge update across warped geometry, symmetry, interaction and surface-ID cases; fixed/moving geometry reuse; repeated near-field buffer clearing; rollback and failed-commit isolation; and cached linear effective systems.
- [compare_coupled.jl](compare_coupled.jl): compare short strong and loose coupling runs against the original source for tip displacement/twist, final structural state, circulation, wake vertices and nodal loads; check exact coupling iteration counts/time/wake rows. Repeat with detailed histories disabled and check the virtual-work identity independently.

Commands from the repository root, using the installed Julia 1.11.5:

```powershell
julia --startup-file=no --compiled-modules=existing --pkgimages=existing --project=. test/efficiencyTests.jl
julia --startup-file=no --compiled-modules=existing --pkgimages=existing --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/test/runtests.jl
julia --startup-file=no --compiled-modules=existing --pkgimages=existing --project=. docs/reviews/efficiency_2026-10-19/compare_coupled.jl compare
julia --startup-file=no --compiled-modules=existing --pkgimages=existing --project=. docs/reviews/efficiency_2026-10-19/measure_implementation.jl
```

The cache flags avoid writes to package caches outside the workspace on this machine. They do not change the package's Julia 1.9 compatibility declaration. Only Julia 1.11.5 was exercised here.

The coupled fixture uses 3 by 6 panels, six physical steps at 30 m/s, 0.01 rad angle of attack, a nonzero 0.05 N pulse, and up to six wake rows. It compares floating-point outputs with `rtol=2e-8, atol=2e-10`; no coupling iteration count changes are accepted. It checks numerical equivalence over a short interval, not long-term LCO amplitude or physical-model validation.

Final verification results:

| Check | Result |
| --- | ---: |
| AeroBeams targeted efficiency, spring and strip-AD regressions | 78 passed |
| Complete WingPropellerUVLM test entry point | 773 passed |
| Coupled original-source equivalence | 18 passed |
| Reduced-history coupled equivalence | 20 passed |
| Independent Pazy virtual-work identity | 1 passed |
| Julia syntax parsing of all changed/new Julia files | 31 files passed |
| `git diff --check` with CRLF handling | Passed |

All 890 executed assertions passed. The repository-wide AeroBeams test suite and full LCO speed sweep were not run. Animation rendering and peak retained process memory were not measured.

## Measurement method

[measure_implementation.jl](measure_implementation.jl) reuses the original diagnostic fixtures/helpers and writes separate `measurements_after.toml` and `cpu_profile_after.txt` files. The original measurement script can be included for helpers and refuses to overwrite the saved baseline when run directly. The coupled capture command also refuses to replace an existing reference.

Each operation is warmed, then sampled seven times using `@timed`. Results record median/minimum time and median allocated bytes, Julia version and thread counts. The UVLM fixture retains 200 bound panels, 20-row wake capacity and the original flight condition. Full wake measurements follow actual accepted steps. The structural fixture uses the same 15-element Pazy model. Active-load checkpoint measurements additionally show the horizon-independent copy behavior for configured loads.

These are bounded allocation diagnostics. Timing is subject to OS scheduling, other running Julia processes, GC and CPU contention. Compare allocation volume primarily; do not interpret a local kernel ratio as a complete simulation speedup. Detailed structural history intentionally retains its required state copies when requested.

## Recorded allocation results

Raw current results: [measurements_after.toml](measurements_after.toml); current trial profile: [cpu_profile_after.txt](cpu_profile_after.txt). Byte counts below are medians per call unless a row describes retained table payload. Decimal MB is used in summaries.

| Operation | Original bytes | Implemented bytes |
| --- | ---: | ---: |
| Inactive distributed-load table payload, 4,096-step Pazy model | 23,598,720 | 0 |
| Newton model checkpoint, 256 steps | 2,311,467 | 811,555 |
| Newton model checkpoint, 4,096 steps | 24,429,771 | 811,555 |
| Scale/interpolate one 4,097-sample load history | 196,967 | 192 |
| Pazy reverse load transfer | 1,202,064 | 0 |
| Pazy forward geometry update | 325,685 | 11,808 |
| AeroBeams assembly, 15 elements | 1,109,648 | 1,023,968 |
| UVLM trial geometry update | 66,012 | 208 |
| UVLM trial, 200 panels and 20 active wake rows | 1,457,950 | 666,405 |
| UVLM near-field force calculation | 289,356 | 267,040 |
| Save requested complete structural step | 91,104 | 92,423 |

The old Newton checkpoint was a full model `deepcopy`; the new row measures the production `checkpoint_model`. The unloaded model's own `deepcopy` now allocates 821,651 B at both horizons after removing empty history payload. A checkpoint of an actively loaded model allocates about 811.5 KB at both horizons as well, while copying that model fully grows from 1.01 MB to 3.77 MB. Appending 4,096 orientation outputs leaves checkpoint allocation about 818 KB. These measurements confirm that the checkpoint no longer copies histories proportional to the simulation horizon.

The 4,096-step checkpoint allocation fell about 96.7%, reusable Pazy geometry allocation about 96.4%, and the full-wake UVLM trial allocation about 54.3%. Reverse load transfer allocates zero in this fixture. `copy_system!` still allocates 30,496 B per reset in this measured backend; it removes most full-System copy allocation, without claiming zero-allocation reset. System storage itself grows from 682,334 B to 687,238 B per System after adding the velocity buffer.

The complete structural saver allocates slightly more because histories now own their state vectors. Response-only Pazy runs avoid invoking it at all. Wake-velocity allocation is essentially unchanged (about 104 KB); the direct induced-velocity arithmetic is retained.

### Timing limits

The timing samples do **not** establish an overall runtime improvement. Other Julia processes were active. The full-wake trial median was approximately 62.0 ms in the original run and 85.9 ms in the current run; near-field force time was 11.5 ms versus 16.4 ms. Even the essentially unchanged wake-velocity operation rose from 47.6 ms to 55.9 ms. These runs do not isolate whether code effects or machine conditions account for the differences. They must not be reported as a full-simulation speedup.

Some isolated routines did show lower times (Pazy load transfer approximately 396 to 9 microseconds; geometry approximately 98 to 23 microseconds), but allocation reductions are the reliable result of this run. A controlled coupled timing comparison, with exclusive resources and the full wake, remains necessary. The current profile continues to place most trial work in induced velocities for near-field forces and normal-velocity wake downwash.

## Follow-up work requiring separate validation

The review's larger proposals remain documented rather than being folded into this patch:

- A typed strip-aerodynamics/ForwardDiff workspace replacing reflective dual-number repair. The current AD path is retained and a selected reference case checks it.
- A complete in-place Newton state restore inventory. This implementation shares only read-only histories and keeps conservative mutable graph copying.
- Preassembled sparse Jacobian structure and direct nonzero-value assembly. The current cache checks the actual pattern safely.
- Cross-surface shared wake-vertex ownership and velocity reuse. Finite-core and interaction policies can differ by receiving surface, so coincident coordinates alone do not justify reuse.
- Parallel induced-velocity receivers, fast summation, matrix-free circulation or GPU kernels. The induced-velocity arithmetic remains the dominant large-wake cost.
- Broader concrete-field changes, cached BC string parsing, general adapter interpolation/topology caching, sparse linear structural assembly, and optional plotting package extensions.
- Application-wide benchmarking on representative long coupled runs, including early/late retained memory, full wake commit cost and LCO amplitude repeatability.

No global BLAS threading policy is imposed. The original 200-panel fixture favored one BLAS thread locally; larger systems and concurrent analyses need their own measurements.

## Julia guidance applied

The implementation follows the official [Julia performance tips](https://docs.julialang.org/en/v1/manual/performance-tips/) on function scope, concrete storage, reducing allocations, views and fused broadcasting, and [LinearAlgebra](https://docs.julialang.org/en/v1.11/stdlib/LinearAlgebra/) guidance on `mul!`, `ldiv!` and sparse factorization reuse. Small fixed vectors are used locally, consistent with [StaticArrays guidance](https://juliaarrays.github.io/StaticArrays.jl/stable/). The broader AD proposal follows [ForwardDiff's configuration/buffer guidance](https://juliadiff.org/ForwardDiff.jl/stable/user/advanced/) and is explicitly deferred. The original review contains the complete source list.
