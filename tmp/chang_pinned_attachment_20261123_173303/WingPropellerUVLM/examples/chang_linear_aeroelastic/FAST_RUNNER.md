# Optional faster Chang runner

This runner uses the existing `chang_case.jl`, configuration loader, structural
integrators, aerodynamic load operator, and time-marching solver. It adds its
own wake-context type and commit method. It does not replace existing methods
or require changes to existing source files.

From the repository root, start a new Julia process with four threads:

```powershell
julia --threads=4 --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/run_chang_linear_aeroelastic_fast.jl
```

The new runner defaults to parallel wake convection and no full per-step
surface geometry history. Wake shedding, vortex stretching, finite-core laws,
interaction groups, symmetry, time step, mesh, and structural/aerodynamic
equations follow the existing case. Parallel work assigns different receiving
vertices to workers; each vertex retains the original order of source
contributions. Near-field forces and the circulation solve remain as before.

Result files are written into a newly allocated directory under
`<configured output directory>/fast_runs/`. Existing output files are preserved.
The runner prints the resulting directory after marching. Animation retains
its separate frame storage when enabled in the case.

Additional optional environment settings:

| Setting | Default | Effect |
|---|---|---|
| `CHANG_SAVE_SURFACE_HISTORY` | `false` | Set `true` to retain every accepted surface geometry. |
| `CHANG_THREADED_WAKE` | `true` | Set `false` to use the original serial wake commit. |

All existing `CHANG_*` case overrides also apply. Setting Julia's thread count
requires starting Julia with `--threads`; a running one-thread REPL cannot
gain worker threads. Existing simulations continue with their loaded code.
The original entry point continues to use its original behavior.

For an already started multithreaded REPL:

```julia
include("lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/run_chang_linear_aeroelastic_fast.jl")
chang_run = run_chang_fast_example()
results = chang_run.results
```

Run the standalone correctness checks with:

```powershell
julia --threads=4 --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/test/chang_runtime_optimizations.jl
```

The checks compare original, optional serial, and threaded wake calculations,
short coupled histories with both structural integrators and coupling schemes,
history retention, animation frames, and output preservation. They also cover
different per-propeller wake capacities, final geometry at a limit stop, and
shared trailing edges with unequal active wake lengths. A direct-induction
check verifies the receiving coordinate beyond a shorter shared-edge wake.
These checks have their own test entry point.

The four-thread check completed with 1,245 assertions passing. A bounded benchmark
on the current 600-panel mesh and 6,760-panel synthetic wake measured median
wake-commit times of 6.810 seconds for the serial code and 2.413 seconds with
four threads, a 2.82x wake speedup. Wake velocities and convected/shed panels
matched exactly in all three repetitions. Other running simulations affect
wall time; this is not an end-to-end speedup measurement. Disabling full surface
history also avoids approximately 2.2 GiB of retained panel data for the current
five-second case.

To repeat the bounded benchmark:

```powershell
julia --threads=4 --project=lib/WingPropellerUVLM lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/validation/benchmark_chang_runtime_optimizations.jl
```
