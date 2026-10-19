# Plot measured/estimated Pazy flutter boundaries without rerunning the sweep.
# Copy the final onset and offset speeds for each pitch angle from
# output/stability_sweep/damping_estimation.xlsx into pazy_boundaries below.

let root = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
    root in LOAD_PATH || pushfirst!(LOAD_PATH, root)
end
using Plots, DelimitedFiles, ColorSchemes
include(joinpath(@__DIR__, "PazyStabilitySweepTools.jl"))

# One row per pitch angle: (angle_deg, onset_speed_mps, offset_speed_mps).
# Use NaN when a boundary has not been identified (or no offset exists).
# Replace the NaN values with your final values from damping_estimation.xlsx.
pazy_boundaries = [
    (0.0, 64.96205412, NaN),
    (1.0, 59.88610453, 75.2735396),
    (2.0, 53.88843598, 65.09014099),
    (3.0, 49.17173599, 57.67919859),
    (4.0, 46.0331531, 53.44946121),
    (5.0, 43.40766244, 50.00245634),
    (6.0, 41.35030159, 46.46354585),
    (7.0, 39.86293642, 44.21992086),
]

pitch_angles_deg = first.(pazy_boundaries)
boundary_rows = Any[]
for (angle, onset_speed, offset_speed) in pazy_boundaries
    isfinite(onset_speed) && push!(boundary_rows, [angle, onset_speed, "onset"])
    isfinite(offset_speed) && push!(boundary_rows, [angle, offset_speed, "offset"])
end

# Styling follows the AeroBeams comparison figure.
aerobeams_color = get(ColorSchemes.darkrainbow, 0.05)
aerobeams_styles = (:solid, :solid)
experimental_colors = (:red, :green)
literature_colors = (:gold, :brown, :magenta)
uvlm_color = :red
uvlm_fill_alpha = 0.25
comparison_labels = (
    uvlm_region="Aerobeams (UVLM)",
    aerobeams_onset="AeroBeams", aerobeams_offset="",
    test_onset_up="Test onset up", test_offset_up="Test offset up",
    test_onset_down="Test onset down", test_offset_down="Test offset down",
    umnast_loss="UM/NAST (Exp. loss)", umnast_panel="UM/NAST (Panel coeffs.)",
    sharpy="Sharpy",
)

ts = 10
fs = 16
lfs = 9
lw = 2
ms = 10
ms2 = 4
msw = 0

comparison_plot = plot_sweep_comparison(boundary_rows, pitch_angles_deg;
    aerobeams_color, aerobeams_styles, literature_colors,
    experimental_colors, uvlm_color, uvlm_fill_alpha, labels=comparison_labels,
    line_width=lw, marker_size=ms, literature_marker_size=ms2,
    marker_stroke_width=msw)
plot!(comparison_plot;
    xlabel="Airspeed [m/s]", ylabel="Root pitch angle [deg]",
    xlims=[30, 90], ylims=[1, 7.25],
    xticks=collect(30:10:90), yticks=collect(0:1:7),
    legend=:topright, tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)

output_directory = joinpath(@__DIR__, "output", "stability_sweep")
mkpath(output_directory)
savefig(comparison_plot, joinpath(output_directory, "pazy_flutter_boundaries.png"))
savefig(comparison_plot, joinpath(output_directory, "pazy_flutter_boundaries.svg"))
display(comparison_plot)
println("Flutter boundary plots saved to $output_directory")
