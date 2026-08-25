# WingPropellerUVLM

`WingPropellerUVLM` is the in-repository Julia package containing the modified
VortexLattice backend and the reusable wing--propeller UVLM routines used by the
AeroBeams coupling work.

The package is intentionally separate from the `AeroBeams` module. This keeps
the global free-wake aerodynamic state independent of AeroBeams' element-local
strip-theory states while allowing both packages to live and evolve in the same
repository.

## Local development

From the repository root:

```julia-repl
pkg> activate lib/WingPropellerUVLM
pkg> instantiate
pkg> test
```

To make the package available in another active Julia environment:

```julia-repl
pkg> develop path="lib/WingPropellerUVLM"
```

Then load it normally:

```julia
using WingPropellerUVLM

grid,ratio = wing_to_grid(xle,yle,zle,chord,theta,phi,ns,nc)
_,ratio,surface = grid_to_surface_panels(grid; ratios=ratio)
system = System([surface]; nw=[maximumWakeRows])
system.reference[] = Reference(Sref,cref,bref,rref,Vref,rho)
snapshot = snapshot_uvlm(system)

advance_uvlm_trial!(system,snapshot,freestream,Δt;
    repeatedPoints=repeated_trailing_edge_points(system.surfaces),
    activeWakeRows=activeWakeRows,
)

# Dimensional body-frame loads for transfer to AeroBeams:
surfaceNodalForces = imperial_nodal_forces(system)
```

The coupled AeroBeams controller should restore the same snapshot before every
outer iteration and call `commit_wake_rows!` only after the structural and
aerodynamic interface residuals converge.

## Near-field force convention

Near-field loads follow Imperial College London's C++ UVLM implementation:
Joukovski segment forces, `-rho*A*normal*gamma_dot` unsteady panel forces, and
the same force-to-vertex transfer including its trailing-edge rule. Fluid
density is stored in `Reference`; the legacy five-argument constructor remains
valid and selects `rho=1.0`.

`imperial_nodal_forces(system)` is the coupling-facing API. It returns one
matrix of dimensional `SVector{3}` body-frame forces per aerodynamic surface.
Freestream force derivatives have not yet been ported to this formulation and
are therefore disabled when Imperial near-field loads are requested; use
`derivatives=false` for the partitioned time-domain coupling.

## Source layout

- `src/backend`: modified VortexLattice implementation;
- `src/wing_propeller`: reusable grid, blade, kinematics, and initialization
  routines from `Wing_Propeller_UVLM`;
- `src/UVLMState.jl`: snapshot, restore, trial, and wake-commit operations;
- `test`: package smoke and transaction tests.

Legacy structural matrices, plotting scripts, and parameter sweeps are not
loaded by this package. The nonlinear structure will be supplied by AeroBeams.

The rigid rectangular-wing free-wake verification case is documented in
[`examples/README.md`](examples/README.md). It uses the same explicit
`propagate_system!` accepted-state loop as the `run_chang*` aeroelastic drivers
and produces ParaView time-series files plus mean spanwise lift results.

See `PROVENANCE.md` and `THIRD_PARTY_NOTICES.md` for source attribution.
