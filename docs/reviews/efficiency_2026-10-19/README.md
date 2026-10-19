**AeroBeams + UVLM efficiency review — 19 October 2026**

The strongest immediate opportunities are reducing copies of full structural load histories, storing only requested results, and making coupling transfers use their two-node interpolation support. Within UVLM, the measured aerodynamic work is concentrated in induced velocity for near-field forces and wake downwash. Dense LU is a smaller part of this particular 200-panel case. Optimizing the existing kernels and eliminating repeated work should precede a new solver architecture.

This review follows the computational paths through `src/`, `lib/WingPropellerUVLM/src/`, the Pazy coupling, the general wing/propeller adapter, and their drivers. It also inspects the strip-aerodynamics/AD path and linear structural integrators, which the current Pazy UVLM example does not use. These paths have different priorities. This is not a validation of every physical model or example in the repository.

Production source files were not changed. The accompanying diagnostic script contains an isolated load-transfer prototype and equivalence checks. The user's existing staged `run_pazy_LCO.jl` was left intact.

**Evidence and reproduction**

Run the bounded diagnostic from the repository root:

```powershell
julia --startup-file=no --project=. docs/reviews/efficiency_2026-10-19/measure_allocations.jl
```

The reviewed machine uses Julia 1.11.5, one Julia thread and six BLAS threads. To avoid writing new package caches outside the workspace, the recorded run additionally used `--compiled-modules=existing --pkgimages=existing`. These particular cache flags are for the installed Julia version; they are not a proposed change to the package's Julia 1.9 compatibility.

`measurements.toml` records seven warmed observations per operation, reporting median elapsed time and allocated bytes. They are bounded diagnostics rather than a controlled performance benchmark suite. Compilation is warmed before measurement, but scheduling, GC, CPU contention and timer resolution still affect the timings. The source contains `@timed` measurements rather than BenchmarkTools estimates. The allocation differences are more reliable than small timing differences.

The aerodynamic fixture has one rectangular symmetric wing, 10 chordwise by 20 spanwise panels, 20 allocated wake rows, 30 m/s freestream and 0.01 rad angle of attack. A full wake is obtained through actual accepted aerodynamic steps. Repeated trials all start from one accepted state. Structural measurements use the actual 15-element Pazy material/geometry data, zero distributed loads and a root clamp. The full coupled LCO sweep was not run, so these numbers are not whole-simulation speedups.

The raw results and the CPU profile of repeated trials with a full wake are in [measurements.toml](measurements.toml) and [cpu_profile.txt](cpu_profile.txt). The runnable diagnostic is [measure_allocations.jl](measure_allocations.jl).

Representative observations, rounded to avoid implying more precision than the sampling supports:

| Operation | Median time, approximately | Allocated per call, approximately |
| --- | ---: | ---: |
| Copy unloaded Pazy model with 256 steps | 2.0 ms | 2.31 MB |
| Copy unloaded Pazy model with 4,096 steps | 9.3 ms | 24.43 MB |
| Assemble AeroBeams system, 15 elements | 1.0 ms | 1.11 MB |
| Scale 4,097 samples of one load history before interpolation | 76 microseconds | 197 KB |
| Interpolate first, then scale the six selected values | 2 microseconds | 608 B |
| Original Pazy load transfer | 396 microseconds | 1.20 MB |
| Preallocated load transfer using nonzero support | 3.3 microseconds | 0 B |
| Original Pazy geometry construction | 98 microseconds | 326 KB |
| Full UVLM trial with 2 active wake rows | 21 ms | 1.46 MB |
| Full UVLM trial with 20 active wake rows | 62 ms | 1.46 MB |
| Wake velocity evaluation with 20 active rows, without convection/shedding | 48 ms | 104 KB |
| Full UVLM System copy, 20-row capacity | 135 microseconds | 700 KB |
| Full AIC assembly, 200 panels | 5.3 ms | 0 B |
| Additional trailing-edge coefficient update | 0.76 ms | 0 B |
| Dense circulation LU/solve, six BLAS threads | 0.46 ms | 322 KB |
| Dense circulation LU/solve, one BLAS thread | 0.35 ms | 322 KB |

The profile is for aerodynamic trials only; it excludes structural solves and wake commit. In the recorded profile, near-field force calculation accounts for roughly two thirds of the trial samples, wake downwash for about one quarter, and AIC assembly for most of the remainder. The exact proportions can change with wake length, mesh, deformation and threading. The much smaller System-copy time changes the priority suggested by a static inspection of `deepcopy` alone.

**1. Highest priority: stop copying future structural load histories for every Newton solve**

Evidence: `src/Problem.jl:659–752`, `src/SystemSolver.jl:110–114`, `src/Model.jl:45`, and `src/Problem.jl:1486–1504`.

`precompute_distributed_loads!` allocates eight arrays of size `(2,3,sizeOfTime)` for every element even when no distributed load of that kind exists. Pazy's aerodynamic loads are nodal BCs, so the arrays remain zero. `solve_NewtonRaphson!` then deep-copies the complete model at the start of every load step, including these arrays. Strong coupling calls that structural solve repeatedly within each physical step. The copied model also contains the growing `R_A_ofTime` history.

For 15 elements, the zero-load arrays occupy `8 categories * 2 nodes * 3 components * 8 bytes * 15 elements = 5760 bytes` per time sample. With the current LCO settings (`duration=10 s`, `U*dt/chord=0.25`, chord `0.0989 m`), that is approximately 66.65 MiB at 30 m/s and 199.96 MiB at 90 m/s. These are calculated storage sizes, not measurements of complete 10-second simulations. The diagnostic directly confirms that the allocated arrays are zero and that model-copy allocation grows with the horizon.

Proposed changes:

- Allocate histories only for load categories that are configured; retain empty arrays of the existing dimensionality for inactive categories.
- Separate immutable model/setup data and precomputed load tables from the mutable state needed to retry Newton. Snapshot the latter into reusable storage. Include the solver's mutated element states, rates, rotations, load caches, constraint state and BC caches in a deliberate restore inventory.
- Keep result history outside that retry snapshot. For a fixed reference frame, use a constant orientation or store orientations only at requested output times, adapting animation indexing accordingly.
- Restore state in place so geometry/load adapters can retain references to stable buffers.

When the simulation horizon has `Nt` samples and average coupling count `K`, copying full `O(E*Nt)` load tables on every step creates an avoidable `O(K*E*Nt^2)` data-copy term, before considering numerical solves. The growing orientation history adds another term of the same temporal order. Replacing the snapshot must preserve the existing retry behavior and cannot be accomplished safely by deleting the copy alone.

**2. Interpolate active distributed loads before scaling them**

Evidence: `src/Core.jl:276–357` and `src/Core.jl:490–517`.

Each active category currently calls expressions such as `interpolate_distributed_loads(problem, σ*f_A)`. Multiplication allocates and scales the entire future history at every element assembly, although the function only needs the current pair of nodal values or two adjacent samples. This affects simulations with applied distributed loads, including AeroBeams cases outside the present Pazy UVLM example.

Instead, pass the load factor into the interpolation routine and scale only its output. A prototype that uses the existing interpolation unchanged, then scales the two three-component outputs, agrees with scaling first for the checked data and changes allocation from about 197 KB to 608 B for a 4,097-sample history. A direct in-place interpolation can remove the remaining small output allocations. Keep the current time/index convention when making this change.

**3. Make saved output proportional to the requested output**

Evidence: `PazyWingUVLMCoupling.jl:99–105,331` and `src/Problem.jl:1606–1675`.

The Pazy driver enables full structural history and explicitly saves it on every accepted step. `save_animation_history=false` disables wing/wake frame capture but does not disable structural snapshots. The saver copies elemental states, nodal states, complementary states, rates and aerodynamic variables; the diagnostic observes about 90 KB per saved step even without strip aerodynamics. At 36,401 samples that suggests roughly 3.3 GB of allocation volume for these snapshot operations alone. Retained memory and allocation volume are different quantities; the full run was not measured.

Add separate controls for scalar response history, detailed structural snapshots, nodal aerodynamic load history and animation frames. Keep tip motion and coupling residuals at every step if needed; make full state snapshots optional or decimated. The LCO driver consumes tip twist and time, so its solver output can be much smaller. Allocate a fixed number of output slots or use `sizehint!` when the count is known.

There is also a correctness issue: `push!(xOverTime,x)` stores the active state vector itself. The diagnostic confirms `xOverTime[end] === problem.x`. While the vector identity remains unchanged, later Newton updates mutate earlier saved entries. `save_load_factor_data!` has the same pattern for `xOverσ`. Save owned copies when preserving these histories, or move to an output layout with explicit independent columns. Avoid replacing necessary snapshot copies with shared mutable references.

**4. Use the two-node coupling support directly**

Evidence: `PazyWingUVLMCoupling.jl:108–117,415–477`.

Each row of `weights` has at most two nonzero entries, but both geometry interpolation and reverse load transfer scan all beam nodes. The loops also create small array slices and temporary cross products for zero-weight nodes.

Precompute the two node indices and interpolation fractions once. Implement `wing_geometry!` and `beam_loads!` using reusable arrays and fixed-size three-component vectors or scalar arithmetic. The reverse mapping must continue to use each beam node's own vortex offset. The same interpolation support must drive the forward and reverse maps so virtual work and moments remain consistent.

The review-only prototype evaluates the exact existing weights, skips their zero entries, reads offsets without slices and fills a caller-owned output matrix. For nontrivial signed forces and offsets, its force/moment output matches the current routine exactly in the recorded fixture. Allocation drops from 1,202,064 B to zero. The isolated kernel is about 120 times faster in this small fixture; the complete coupled solver will improve by a much smaller amount because aerodynamic and structural solves remain.

This change is also simpler: two indices and two weights describe the actual mapping more directly than a mostly-zero dense matrix plus loops across every beam node.

**5. Remove repeated aerodynamic work before replacing the algorithm**

Evidence: `backend/analyses.jl:735–810`, `backend/induced.jl:430–705,823–1027`.

The dynamic backend unconditionally sets `surface_motion=true`, assembles the complete AIC using the current wake-shedding locations, and immediately recomputes the trailing-edge columns using the same geometry/options. The diagnostic found an exact zero change from the second update in four single-surface configurations combining symmetry and finite-core flags. The full assembler's source already includes those trailing-edge/interface contributions.

Proposed control flow:

```julia
if full_matrix_required
    influence_coefficients!(AIC, surfaces; coefficient_options...)
elseif shedding_geometry_changed
    update_trailing_edge_coefficients!(AIC, surfaces; coefficient_options...)
end
```

Before changing production code, compare deformed and multiple-surface cases, interaction masks, trailing-vortex settings and custom core rules. Preserve the trailing-edge-only path when the full matrix was not rebuilt. A stationary surface does not by itself prove that AIC or its factorization is reusable: shedding locations and relevant solver options also matter.

The full AIC kernel already has negligible heap allocation in the fixture. Its cost is arithmetic and memory access, so scattering additional `!` annotations over it would not address the main issue. Its receiving-panel loop writes across matrix rows, trading column-major output access against edge reuse and source traversal; any blocked layout change should be measured rather than applied mechanically.

The current coupling performs an additional UVLM evaluation when final geometry first meets the convergence tests. That evaluation checks load consistency at the geometry being accepted. Retain this invariant. Reducing it requires a revised acceptance algorithm and validation, rather than simply omitting the call.

**6. Concentrate UVLM kernel work on near-field forces, downwash and wake velocity**

Evidence: `backend/nearfield.jl:77–133,177–263`, `backend/circulation.jl:119–188`, `backend/wake.jl:302–477` and `cpu_profile.txt`.

Near-field forces and wake downwash dominate the recorded full-wake trial. Wake velocity evaluation is separately expensive during commit; the coupling already advances the wake once per accepted physical step, which should be preserved.

Allocate the near-field force/velocity buffers once and overwrite them. The default `near_field_forces!` currently creates new chord, span, unsteady and velocity matrices, plus a negated circulation array, on each call. Although the System already owns force buffers, `propagate_system!` only supplies those as reusable inputs to the legacy force functions. Add a uniform workspace interface to the active implementation. This reduces allocation; the induced-velocity arithmetic will still require most of its time.

Cache immutable topology: surface circulation offsets, interaction sender lists and known shared trailing-edge connectivity. `repeated_trailing_edge_points` currently performs a pairwise geometric search when invoked, including before a trial that defers wake advancement. It should be computed where needed, with a valid cache rule. Precomputing arbitrary geometric coincidences is unsafe when independent surfaces can move into or out of coincidence.

Parallel candidates are receiving control points, force segments, and wake vertices, since each can write to an independent destination. AIC edge reuse writes to several source columns in one receiving row: parallelize receiving rows or disjoint receiving blocks, keeping the source-edge accumulation within a row local. Wake duplicate-point reuse needs explicit owner/copy ordering before threading.

In the active `get_wake_velocities!`, the `continue` after copying a duplicate point is inside the loop over duplicate matches. It does not skip the enclosing receiving-vertex computation, so the intended reuse is lost. This is a static control-flow finding for shared-edge, multiple-surface cases; it is not exercised by the single-wing benchmark. Resolve ownership and skip the correct outer computation, checking shared-point velocities before and after.

For sufficiently large meshes, fast summation or a matrix-free iterative formulation can be worthwhile. The direct interaction costs scale approximately with `Nb^2`, `Nb*Nw` and `Nw*(Nb+Nw)` for the relevant stages, while dense factorization scales as `Nb^3`. Such methods add numerical approximation, convergence and maintenance costs. The current 200-panel measurements favor the changes above first.

**7. Reuse structural factorization at the correct lifetime**

Evidence: `src/SystemSolver.jl:148,285–324`, `src/Core.jl:1181–1286`, `src/Problem.jl:611–619`, and `aeroelastic/StructuralTimeIntegration.jl:185–202,357`.

AeroBeams already uses a sparse Jacobian. The issue is repeated indexed block insertion and repeated factorization, rather than a missing sparse solver. `skipJacobianUpdate` skips Jacobian calculation, but `solve_linear_system!` still evaluates a fresh backslash solve. Cache the numeric factorization while the matrix is unchanged. When values change but the structural pattern is fixed, reuse the sparse symbolic factorization. Precompute a conservative sparsity pattern from connectivity and DOF layout; discovering the pattern only from a zero-deformation numerical matrix can miss entries that become nonzero later.

The expression `-jacobian\residual` negates the matrix before solving. Factor the unnegated Jacobian and solve a negated RHS into `Δx`. For the existing `ldiv!` API, provide a separate RHS buffer if required by the chosen factorization. Invalidate numeric caches after Jacobian rebuilds, retries, topology/constraint changes and problem replacement. Preserve the current assignment/accumulation rules for special nodes.

The newer linear `partitioned_aeroelastic_step` already factors its effective matrix once per physical step and reuses it within the coupling iterations. This is good. If `M,C,K,dt` and method parameters stay constant across the run, move that cache to run scope. If they change, rebuild it. The older `partitioned_generalized_alpha_step` rebuilds its effective system through each corrector and duplicates the newer iteration loop; a compatibility adapter to the common loop would remove both divergence and repeated work after acceptance/residual semantics are compared.

**8. Reduce small-array churn in AeroBeams element kernels**

Evidence: `src/Core.jl:134–224,229–357,986–1180` and `src/Utilities.jl:65,162`.

Examples include `u .= x[DOF_u]`, matrix products on the RHS of `.=` assignments, `Matrix(RR0')`, `Matrix(R_p1')`, `mul3` building three products then concatenating them, and fresh zero arrays for absent loads. `.=` controls the destination; it does not make ordinary matrix products or indexed RHS slices allocation-free.

Use views or explicit short gathers for state extraction, fused broadcasting for elementwise formulas, `mul!` with owned scratch for dynamic matrices, and transpose wrappers where consumers support them. For fixed three-vectors and 3-by-3 matrices, selective `SVector`/`SMatrix` kernels are natural. The diagnostic observes about 35 KB allocated by one rotation-variable update and 1.1 MB by one 15-element system assembly at zero deformation. These routines are repeated within Newton and coupling.

Some rotation derivatives are calculated before `skipJacobianUpdate` causes the caller to return. Separate rotation quantities required for residuals/rates from derivatives required only for Jacobian updates. Trace dependencies before skipping a derivative, since transient rates also use some rotation derivatives.

Use static arrays selectively: globally encoding mesh sizes or large matrices in types can increase compilation and code size. Keep runtime-sized state vectors, AIC matrices and large beam matrices in normal or sparse arrays.

**9. Repair type uncertainty at numerical boundaries**

Evidence: `src/Element.jl:363–367,435–451`, `src/AeroProperties.jl:18–34,66–74`, `src/AeroSurface.jl:10–27`, `src/Problem.jl:376`, and `backend/system.jl:93`.

Several hot structures store abstract fields (`Real`, `Function`, arrays with unspecified dimensionality), and strip-flow structures have unannotated fields, hence `Any`. In contrast, accepting an abstract argument type in a function signature is usually compatible with specialization; it is the stored field or container representation that needs attention.

For structural geometry that is always Float64 in the current implementation, normalize values once into concrete storage. Declare known load-history dimensionality, for example `Array{Float64,3}`. Isolate numeric kernels behind function boundaries so the compiler sees concrete values and callables without parameterizing the complete cyclic Beam/Element/Model graph.

The strip-flow state deliberately receives ForwardDiff dual numbers. Replacing its fields with Float64 would break AD. Use typed primal and AD workspaces, or a typed local aerodynamic kernel, so differentiation does not mutate a large graph and then require reflective repair. Diagnose actual specialization with `@code_warntype` and allocation profiles before broad type redesign.

**10. Simplify the AD and boundary-condition update paths**

Evidence: `src/Core.jl:844–981,2131–2160`, `src/BC.jl:104–142`, and `PazyWingUVLMCoupling.jl:370–380`.

The strip aerodynamic derivative routine concatenates state arrays, allocates a Jacobian, mutates aerodynamic structures with dual values, then recursively visits fields to reset them. Prefer a numerically typed kernel with owned outputs and reusable ForwardDiff configuration/Jacobian buffers. This applies to strip models and hinge constraints; the Pazy UVLM beam has `aeroSurface=nothing`, so it is not a current Pazy priority.

BC updates allocate many six-component arrays and rescan fixed string type descriptions during every load exchange. Parse the BC layout once during construction. A dedicated numeric load update can fill the appropriate cached values directly while retaining `update_BC_data!` as the general path for time-dependent, transformed and follower inputs. The current Pazy code correctly updates BCs on the active model after a Newton retry; retain that behavior.

**11. Simplify setup and optional components without overstating their runtime effect**

| Area | Concrete simplification | Relevance |
| --- | --- | --- |
| UVLM transaction state | Consolidate full-System trial copies and the separate snapshot/restore API around one explicit restore contract; reuse accepted/trial/commit buffers and swap only after successful convection | Repeated allocation and maintenance; measured copy time is small for this fixture |
| Trial geometry | Reuse `update_surface_panels!` with stored ratios instead of allocating via `grid_to_surface_panels` then copying its output | Every aerodynamic trial; retain matching core/normal behavior |
| Structural propeller attachment | Scatter local blocks directly to `node_dofs`, replacing `zeros(6,wing_dof_count)` selectors and dense products through them | Simpler setup and better scaling |
| Linear beam assembly | Use sparse triplets or a band/block representation when mesh size warrants it; avoid retaining redundant expanded matrices unless requested for diagnostics | Alternative linear structural model; one-time setup |
| Spanwise property inputs | Validate and normalize constant/vector/tabulated inputs once, then cache quadrature data and breakpoints | Setup; existing general input support is useful |
| `geometry.jl:137` interpolation | Use component views/direct writes instead of repeated slices and `vcat`; make the reused-y-coordinate assumption explicit | Geometry preprocessing, usually below solver priorities |
| General UVLM adapter | Cache span stations, undeformed chord/xLE distributions and interpolation support; reconstruct full beam kinematics once per trial and reuse it for hubs | Repeated trial setup in wing/propeller cases |
| Quadrature helper | Store fixed Gauss nodes/weights in constant tuples instead of rebuilding a Dict and all rules per call | Strip loads with tip correction |
| Scaling utilities | Rename `divide_inplace!`/`multiply_inplace!`, which return generators of newly scaled values, or replace uses with direct expressions | Clarity; the existing names imply mutation that does not occur |
| Duplicate coupling loops | Retain integrator-specific formulas, share the iteration/acceptance orchestration after checking semantics | Reduces maintenance and inconsistent optimizations |
| Commented implementations | Remove obsolete commented function bodies from circulation, wake, induced and analysis files; use Git history | Readability; no execution-time benefit |
| Legacy force implementation | Keep it only for a documented comparison/compatibility purpose, with a common buffer contract | Distinct physical formulations need comparison before removal |
| Package startup | Move development-only Revise usage out of `AeroBeams.jl`; make visualization loading optional; put LCO orchestration in `main(config)` | Startup, repeated workers and maintainability |

The existing `UVLMSnapshot` is a useful starting point, not a complete replacement for the facade's deep-copy semantics. It omits some System fields, while trial application changes geometry and other caches. Reusing buffers requires auditing every mutated field. A commit failure must continue to leave the accepted state unchanged and the successful trial retryable.

The eager `repeated_trailing_edge_points` calculation in trial evaluation is unnecessary when the deferred wake path does not use it. Sharing immutable topology is a useful simplification. A new generic cache framework for every helper is unlikely to pay for itself.

**12. Broader modes and case-level parallelism**

`src/Problem.jl:798–858` computes the complete dense eigendecomposition and filters to `nModes` afterward. For large stability studies requesting a few modes, investigate a generalized sparse/shift-invert approach that handles the descriptor equations and infinite modes correctly. The current determinant check can also add a factorization; inspect factorization success or appropriate conditioning diagnostics instead. This is a separate scalability project from Pazy time marching. In trim solves, avoid forming an explicit pseudoinverse or an inverse-like map when a rank-aware least-squares factorization can express the intended solution, preserving the current regularization and rank-deficient behavior.

The LCO driver starts independent simulations for each speed. After reducing per-case memory, separate worker processes can run independent speeds concurrently, each with local state and controlled BLAS threading. Keep plotting outside concurrent solver work. The process startup and copied memory costs should be included in any comparison. Increasing Julia's thread count alone will not parallelize these serial loops.

For a 200-by-200 circulation matrix, the recorded one-thread BLAS solve is faster than the six-thread solve; this is a local observation, not a universal setting for larger matrices. Start by comparing one BLAS thread for small cases and when parallelizing cases. The LCO driver already disables UVLM/animation histories, but structural history still needs the separate control described above. Frequent Newton and progress printing in `run_pazy_wing_uvlm.jl` should be decimated for timed runs.

**Suggested implementation order and verification**

1. Skip inactive distributed-load histories and stop scaling full histories inside residual assembly. Confirm no-load, active constant/time-varying and retry behavior.
2. Replace full-model Newton checkpoints with owned state buffers, preserving retries and load/BC state. Separate detailed output from scalar response histories; fix history aliasing.
3. Implement the two-node forward/reverse maps with reusable buffers. Check force, moment and virtual-work equivalence, including deformed geometry and clamped-root symmetry.
4. Remove the redundant trailing-edge pass after a full rebuild, after extending the checked flag/geometry combinations. Make geometry and near-field updates reuse existing buffers.
5. Reuse Jacobian factorization and sparsity metadata, then tune the small-array structural kernels and AD paths as their profiles warrant.
6. Consider parallel receivers/cases and approximate aerodynamic summation only with new profiling evidence.

Use the current mesh, wake, timestep, pulse, coupling tolerances and Newton tolerances for comparisons. Check tip displacement/twist, aerodynamic resultants, convergence residuals, coupling iterations, wake/circulation commit behavior and repeatability. A shorter wake, larger timestep, looser tolerances or loose coupling changes the numerical approximation and cannot establish a code-only optimization gain.

For a complete performance comparison, warm the relevant paths and profile a coupled interval with a full wake, both early and late in a longer run. Record elapsed time, allocation volume, peak retained memory and iteration counts. Preserve raw reference output and compare each change individually. BenchmarkTools setup/evals controls are useful for mutating kernels; Julia's CPU and allocation profilers locate the expensive call stacks. Do not add a broad framework until these small, concrete changes have been measured.

**Official Julia and package guidance consulted**

- [Julia performance tips](https://docs.julialang.org/en/v1/manual/performance-tips/): concrete fields, specialization boundaries, allocation, array access and threading. Apply these selectively; views can be slower for repeatedly used irregular data.
- [Julia LinearAlgebra](https://docs.julialang.org/en/v1.11/stdlib/LinearAlgebra/): `mul!`, `ldiv!` and `lu!(F,A; reuse_symbolic=true)` for fixed sparse patterns.
- [Julia profiling](https://docs.julialang.org/en/v1.11/manual/profile/): CPU sampling and `Profile.Allocs` for allocation origins. Keep compilation out of the measured interval.
- [BenchmarkTools manual](https://juliaci.github.io/BenchmarkTools.jl/stable/manual/): interpolate benchmark inputs, control sample setup, and use one evaluation per sample when mutation requires resetting input.
- [StaticArrays documentation](https://juliaarrays.github.io/StaticArrays.jl/stable/): fixed small vectors/matrices can improve local arithmetic and storage; large static arrays can increase compilation and code size.
- [ForwardDiff advanced usage](https://juliadiff.org/ForwardDiff.jl/stable/user/advanced/): reuse differentiation configuration and owned buffers, choose chunking with measurements, and preserve numeric types needed for dual propagation.
