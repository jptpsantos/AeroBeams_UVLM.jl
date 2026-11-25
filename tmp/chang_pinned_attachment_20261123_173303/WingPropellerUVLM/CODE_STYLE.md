# AeroBeams notation and implementation conventions

The custom code uses the local AeroBeams implementation as its reference,
especially `src/Beam.jl`, `src/Model.jl`, `src/Element.jl`, `src/Core.jl`, and
`src/SystemSolver.jl`. AeroBeams files themselves are excluded from this pass.

The 2026-11-10 pass covers all 211 existing custom Julia files: the UVLM
backend and facade, geometry and kinematics, structural integrators, study
package, general wing–propeller examples, Pazy coupling, both Chang models,
validation scripts, tests, and archived custom scripts. Scratch files and
generated simulation results are not production source.

| Item | Convention |
|:--|:--|
| Types/modules | `UpperCamelCase`, such as `UVLMDynamicProblem` |
| Functions | `snake_case`; existing `create_TypeName` constructors follow AeroBeams |
| Mutation | Retain `!` on established mutating operations |
| Internal descriptive names | `camelCase`, such as `acceptedLoad`, `beamNodes`, `nodeIndex` |
| Main objects | Explicit `problem`, `model`, `workspace`, `aeroSolver`, `operatingPoint` |
| Constants | Established uppercase names and mathematical constants |
| Layout | Four spaces, explicit loops, spaces around binary operators, compact comma-separated lists, `keyword=value` |
| Algorithms | Preserve numerical expressions, accumulation order, tolerances, units and discretization |

Mathematical definitions take priority over a mechanical naming rule:

| Notation | Meaning and constraint |
|:--|:--|
| `r`, `r_ref_A`, `r_spun_A` | Position vectors; the suffix identifies the stated basis/configuration |
| `u`, `udot`, `uddot` | Translational displacement and its time derivatives where those are the actual quantities |
| `p`, `pdot`, `pddot` | AeroBeams Wiener–Milenković rotation parameters and parameter derivatives; `pdot` is not automatically spatial angular velocity |
| `θ`, `ψ`, existing Euler-angle names | Actual Euler/modal angles; retain the implemented rotation order and signs |
| `q`, existing generalized-state names | Mixed structural coordinates; do not rename an entire mixed vector as a translational `u` or a rotation `p` |
| `F`, `M`, `V`, `Ω` | Forces, moments, velocities and angular velocities in their documented contexts |
| Structural `M`, `C`, `K` | Existing assembled linear mass/damping/stiffness matrices; distinct from a moment vector `M` |
| AeroBeams beam `S`, `I` | Sectional stiffness and inertia, passed through the original AeroBeams definitions |
| `Γ`, `dΓdt` | Circulation and its time derivative |
| `_A`, `_b`, `_S`, `_U`, `_n`, `_np1`, `_n1`, `_n2` | Meaningful basis, time-level and node subscripts; retain them |

The UVLM global basis is downstream/span/up. The linear Chang structural
basis is span/downstream/down. Pazy's AeroBeams basis uses its existing
`A_TO_UVLM` conversion. The nonlinear Chang adapter places its AeroBeams model
in the UVLM global basis and separately maps sectional directors. A legacy
UVLM `_A` suffix does not identify AeroBeams basis A merely by sharing a letter.
The interface matrices and physical assumptions are unchanged.

Existing public names are preserved by the user's explicit choice. Therefore
backend properties, constructor keywords, case settings, named-tuple keys,
environment variables, serialized names and CSV columns can retain
`snake_case`. Locals used directly by property/keyword shorthand or generated
diagnostic code may retain the same spelling. Those compatibility boundaries
do not need a translation layer or duplicate API.

The pass standardized 956 descriptive internal names and clarified the object
names in the core facade and study API. Propeller position vectors now use
`r_*` instead of `p_*`. Julia syntax-tree comparisons check the numerical
implementations against the pre-pass snapshot, allowing the recorded
identifier changes and ignoring source-line metadata. The separately reviewed
entry guards and diagnostic/test corrections are listed in the validation
report. Public keys, case values and numerical expressions are preserved in
the implementations. AeroBeams source/test/development file hashes are checked
independently.

Validation results and the source inventory are recorded in
`NOTATION_VALIDATION.md`. Start a fresh Julia process after changing source;
existing simulations keep their already loaded methods, and study provenance
checks still require the source corresponding to their saved metadata.
