# Root-strain spectrograms from saved sweeps

From the repository root, generate the completed 5-degree **upsweep**:

```powershell
julia --project=. lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/plot_pazy_LCO_spectrogram.jl
```

The default input is the newest numbered run inside
`output/LCO_downsweep_pulse/AoA5`, relative to the example directory. Before
numbered runs exist, saved data directly in `AoA5`, then legacy 5-degree data
directly in `output/LCO_downsweep_pulse`, are accepted. The saved settings identify
the plotted incidence; the current simulation driver is not used to infer what
was run. This is postprocessing only, and does not load an aeroelastic solver.
PNG, SVG and processing settings are written into the selected run's
`spectrograms` folder, keeping angles and simulation runs separate.

To process both independently saved branches:

```powershell
julia --project=. lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/plot_pazy_LCO_spectrogram.jl lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/output/LCO_downsweep_pulse both
```

From a Julia REPL or the IDE, include the file and call the function explicitly:

```julia
include("lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/plot_pazy_LCO_spectrogram.jl")
results = plot_pazy_LCO_spectrogram(;
    angle_of_attack_deg=3.0,      # newest saved 3-degree run; use 5.0 or 7.0 as needed
    plot_branches=:upsweep,       # :upsweep, :downsweep, or :both
    window_seconds=1.0,
    overlap=0.95,
    frequency_limits=(0.0, 100.0),
    color_limits=(-40.0, 60.0),
    airspeed_limits=(0.0, 52.0), # optional blank area below simulated speeds
    detrend=:linear,
    reference_strain_ue=1.0,
    save_spectra=true)            # optional long-format spectrum CSV
spectrum = results[:upsweep].spectrum
```

Including the file only defines functions. The bottom of the file contains
editable defaults for command-line execution. Functions accept an explicit
`input_directory` and a separate `output_directory` if desired. An explicit input
directory takes priority over angle-based selection and can select an older run.
For example, pass `.../LCO_downsweep_pulse/AoA7/run_001` as the CLI directory
argument, followed by `both`. Set `angle_of_attack_deg` near the bottom of the
script when using automatic selection from the command line. Missing 3/7-degree
data raises an error rather than selecting the legacy 5-degree histories.

## Relation to Figure 4

The supplied *Effects of Separated-Flow Nonlinearities on the Flutter Behavior
of the Pazy Wing* (dos Santos, Cesnik and Marques, SciTech 2024), Figure 4,
uses **upsweep** root longitudinal strain histories. Panels (d), (e) and (f)
are the 5-degree experimental, nonlinear-model and linear-model results,
respectively. The text describes a primary band around 28 Hz, its integer
harmonics, and a decrease of the dominant frequencies as speed increases.
These are comparison observations, not imposed frequencies in this code.
The saved AeroBeams/UVLM response determines the plotted bands.

The paper does not give the STFT window length, overlap, detrending algorithm,
FFT normalization or dB reference. Consequently this script implements a
documented amplitude convention; it does not claim identical numerical color
levels to Figure 4. AeroBeams/UVLM is also a different model from the paper's
UM/NAST strip-aerodynamics formulations.

## Processing convention

- Read the signed `root_longitudinal_strain_ue` signal, already in microstrain.
  Taking its absolute value before the FFT would change the harmonic content.
- Use uniformly sampled times, a symmetric Hann window of approximately 1 s,
  and 95% overlap. Only complete windows are included; there is no time padding
  or extrapolation at the ends of the sweep.
- Remove a least-squares straight line within each window before tapering.
  This suppresses the changing static bending contribution. It is not an
  independently computed static-equilibrium subtraction, and can affect low
  frequencies. Use `detrend=:mean` or `:none` to inspect that sensitivity.
- Compute the one-sided peak-amplitude spectrum as
  `A = 2*abs(rfft(window .* detrended_strain))/sum(window)`.
  DC and an even-length FFT's Nyquist bin use a factor of 1 instead of 2.
- Plot `20*log10(A/reference_strain_ue)`. By default, 0 dB means a spectral
  amplitude of 1 microstrain, 20 dB means 10 microstrain, and 40 dB means
  100 microstrain. This is an amplitude spectrum, not a power spectral density.
  There is no normalization to each frame's largest peak.
- Pad the FFT to the next power of two at least four times the window length.
  This refines the frequency grid; the physical resolution remains governed
  by the approximately 1 s window (a scale of approximately 1 Hz, with a wider
  Hann main lobe). Longer windows improve frequency discrimination and smear
  changes over a wider range of speeds.
- Associate each spectrum with the airspeed at its window-center time. For
  the saved 40 s upsweep, a 1 s window spans approximately 0.425 m/s; for the
  20 s downsweep it spans approximately 0.85 m/s. Descending frames are sorted
  by speed only after the Fourier transforms, keeping their columns aligned.
- The dB representation is floored at -80 dB; default color limits are
  -40 to 60 dB for both branches. Color clipping does not change stored
  amplitudes. `save_spectra=true` exports time, speed, frequency, amplitude and
  floored dB values; `results[branch].spectrum` exposes the same arrays.

The simulated gauge is the existing first-element midpoint approximation to
root strain. The downsweep includes its applied tip pulse and resulting
transient; the nonzero force samples are recorded in the settings and their
start is marked on the figure. Its broadband content should be interpreted
with that forcing history in mind. The downsweep is an additional diagnostic;
Figure 4 itself uses upsweeps.

Run the calibration and CSV/plot checks without a simulation:

```powershell
julia --project=. lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/test_pazy_LCO_spectrogram.jl
```
