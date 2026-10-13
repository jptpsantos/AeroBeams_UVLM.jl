`moving_block_2.jl` is an unchanged copy of the supplied moving-block script.
Its final lines are a standalone example that expects `t` and `t_prop` to
exist, so the file is kept as a reference rather than included directly.

The callable `moving_block_damping` function is in `../src/MovingBlockDamping.jl`.
It has the original five settings, displays the two diagnostic plots, and
returns growth rate, frequency and angular frequency. `moving_block_metrics`
remains available there for study and sweep diagnostics. The Pazy single-case
runner passes its wingtip twist history as `result.time` and
`result.tip_twist_degrees`.
