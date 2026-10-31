# Interpolated Pazy bifurcation diagrams

Four saved UVLM runs, all at 5 degrees root incidence, plotted separately because their simulation settings differ. Each folder contains a four-panel figure, individual response panels, interpolated CSV data and source notes. Existing sources and saved results were left unchanged.

| Saved run | Four-panel PNG | Vector PDF | Interpolated data |
| --- | --- | --- | --- |
| 20261024_165844546 | [Figure](20261024_165844546/pazy_bifurcation_interpolated.png) | [PDF](20261024_165844546/pazy_bifurcation_interpolated.pdf) | [CSV](20261024_165844546/interpolated_curves.csv) |
| 20261026_111045715 | [Figure](20261026_111045715/pazy_bifurcation_interpolated.png) | [PDF](20261026_111045715/pazy_bifurcation_interpolated.pdf) | [CSV](20261026_111045715/interpolated_curves.csv) |
| 20261026_203829683 | [Figure](20261026_203829683/pazy_bifurcation_interpolated.png) | [PDF](20261026_203829683/pazy_bifurcation_interpolated.pdf) | [CSV](20261026_203829683/interpolated_curves.csv) |
| 20261027_104048381 (latest) | [Figure](20261027_104048381/pazy_bifurcation_interpolated.png) | [PDF](20261027_104048381/pazy_bifurcation_interpolated.pdf) | [CSV](20261027_104048381/interpolated_curves.csv) |

Solid curves interpolate the final-window minimum and maximum envelopes using PCHIP on the envelope center and nonnegative half-range. The curves pass through the saved envelope values without crossing or extrapolating beyond the sampled speeds. Triangles show recorded extrema; black crosses mark envelope samples whose convergence was not verified. The oldest run has no local-extrema CSV, so its triangles show window extrema only. These are envelope interpolations, not continuation or Hopf-point calculations.

Colors follow `dev/sweptPazy/S0/PazyWingPitchRangeURangeLCOPlotGenerator.jl`: `cgrad(:rainbow, 5; categorical=true)` for incidences `[1,3,5,7,10]` degrees. Consequently all saved 5-degree runs are green. Line width is 2.8 pt, steady-curve width is 1.8 pt, and triangle marker size is 3.4 pt. The fourth panel is geometric incidence (root incidence + elastic twist).

To regenerate all saved cases from the repository root, with outputs in a fresh directory:

```powershell
julia --project=. lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/plot_pazy_bifurcation_interpolated.jl --all
```

Omit `--all` for the latest saved run. Pass one summary CSV path to select a run, or several paths to overlay cases. Including the script also exposes `PazyInterpolatedBifurcation.plot_interpolated_bifurcation` with `linewidth`, `markersize`, `angle_of_attack_deg`, axis limits and other appearance options. No simulations are started.
