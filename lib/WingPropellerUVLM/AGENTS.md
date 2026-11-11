# AeroBeams conventions for custom code

Use the repository's AeroBeams `src/Beam.jl`, `src/Model.jl`,
`src/Element.jl`, `src/Core.jl`, and `src/SystemSolver.jl` as the primary
notation and implementation references. Leave the AeroBeams files unchanged.
These instructions cover this package, its studies, examples, and tests.

- Use `snake_case` for functions and `create_TypeName` for new public type
  constructors when applicable. Use `camelCase` for new descriptive internal
  variables. Preserve existing public function names, keyword arguments,
  properties, record keys, environment variables, and output-file schemas.
- Retain mathematical notation and basis/time/node subscripts. Use `r` for
  position vectors and `p` for AeroBeams rotation parameters. Keep Euler angles
  distinct from Wiener–Milenković parameters. Document coordinate conversions
  at the interfaces; preserve their signs, axes, and physical definitions.
- Use four-space indentation, explicit numerical operations and loops, compact
  argument/index lists, and `keyword=value` as in AeroBeams. Keep the established
  module and file organization. Avoid extra wrappers, configuration layers,
  dependencies, or model changes for stylistic uniformity.
- Prefer a local readable correction. Preserve operation order, solver
  tolerances, units, discretization and numerical algorithms during notation
  changes. Run relevant existing regression tests; do not change numerical
  expectations to accommodate a notation edit.

See [CODE_STYLE.md](CODE_STYLE.md) for definitions and compatibility boundaries.
