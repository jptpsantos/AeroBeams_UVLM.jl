# Source provenance

The aerodynamic source in this package was imported from
`jptpsantos/Wing_Propeller_UVLM` at commit
`12f58e73fedc09c32c6fc4a2305fdc62f5a48872` (2026-08-24).

The files in `src/backend` are a modified copy of the MIT-licensed
`byuflowlab/VortexLattice.jl` implementation. They include the wing--propeller
interaction mask and time-domain segmented/unsteady force storage used by the
research model. The original VortexLattice copyright and license are retained
in `LICENSE` and `THIRD_PARTY_NOTICES.md`.

Files in `src/wing_propeller` originate from the wing--propeller research
repository and are now loaded through the `WingPropellerUVLM` module rather than
through top-level `include` calls.
