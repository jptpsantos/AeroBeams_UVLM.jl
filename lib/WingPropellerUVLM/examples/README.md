# Examples

The [`chang_linear_aeroelastic`](chang_linear_aeroelastic/README.md) directory
contains the copied Chang linear structural model coupled to the current
Imperial-load, free-wake `WingPropellerUVLM` implementation.

## Rigid rectangular wing with a free wake

`rectangular_wing_free_wake.jl` follows the explicit accepted-state stepping
pattern used by the `run_chang*` aeroelastic drivers. Each physical step:

1. copies the accepted surface into `previous_surfaces`;
2. snapshots the accepted aerodynamic history;
3. calls `advance_uvlm_trial!`, which restores the snapshot and invokes
   `propagate_system!`;
4. commits exactly one new wake row; and
5. records total and sectional lift coefficients.

The rigid example performs one aerodynamic trial per step. In a partitioned
AeroBeams analysis, steps 2--3 may be repeated for each structural iterate, but
step 4 must occur only after the coupled solution is accepted.

From the repository root, run:

```powershell
julia lib/WingPropellerUVLM/examples/rectangular_wing_free_wake.jl
```

The example activates the local `lib/WingPropellerUVLM/Project.toml`
automatically, so the same file can also be launched with VS Code's **Run Julia
File** command. On a fresh clone, instantiate its dependencies once with
`julia --project=lib/WingPropellerUVLM -e "using Pkg; Pkg.instantiate()"`.

The defaults define a 7.5 m by 1 m rectangular wing at 5 degrees angle of
attack and 100 m/s. The full-span lattice has 30 spanwise by 4 chordwise panels.
The time step is one quarter of the chord-passage time, the simulation lasts 15
chord-passage times, and the reported mean uses the final five chord-passage
times. Air density is 1.225 kg/m^3 and is supplied directly to `Reference` for
the Imperial near-field force evaluation and coefficient normalization.

Outputs are written under `examples/output/rectangular_wing_free_wake`:

- a ParaView `.pvd` time collection and time-indexed `.vtm` files;
- `lift_history.csv`;
- `mean_lift_distribution.csv`;
- `mean_lift_distribution.svg`; and
- `summary.txt`.

With the Imperial near-field force implementation, the default case gives the
following regression values when averaged over the final five chord-passage
times:

| Quantity | Value |
| --- | ---: |
| Mean `CL` from `body_forces` | `0.3894275619460095` |
| Mean `CL` from sectional integration | `0.3894275619460094` |
| Mean lift | `17889.328626894807 N` |

The agreement of the two independently assembled lift coefficients is also
checked by the package test suite at a smaller lattice size.

The following environment variables support convergence and shorter smoke
studies:

| Variable | Default | Meaning |
| --- | ---: | --- |
| `UVLM_RECT_SPAN` | `7.5` | Full span in metres |
| `UVLM_RECT_CHORD` | `1.0` | Chord in metres |
| `UVLM_RECT_ALPHA_DEG` | `5.0` | Angle of attack in degrees |
| `UVLM_RECT_VINF` | `100.0` | Freestream speed in m/s |
| `UVLM_RECT_RHO` | `1.225` | Air density in kg/m^3 |
| `UVLM_RECT_NS` | `30` | Spanwise panel count |
| `UVLM_RECT_NC` | `4` | Chordwise panel count |
| `UVLM_RECT_STEPS_PER_CHORD` | `4` | Steps per chord-passage time |
| `UVLM_RECT_CONVECTIVE_TIMES` | `15` | Total chord-passage times |
| `UVLM_RECT_AVERAGING_TIMES` | `5` | Final chord-passage times averaged |
| `UVLM_RECT_VTK_STRIDE` | `5` | Physical steps between VTK frames |
| `UVLM_RECT_OUTPUT_DIR` | package output directory | Output directory |

For example, a fast smoke run is:

```powershell
$env:UVLM_RECT_NS = "8"
$env:UVLM_RECT_NC = "2"
$env:UVLM_RECT_CONVECTIVE_TIMES = "2"
$env:UVLM_RECT_AVERAGING_TIMES = "1"
julia lib/WingPropellerUVLM/examples/rectangular_wing_free_wake.jl
```
