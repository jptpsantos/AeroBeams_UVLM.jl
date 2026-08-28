const EXAMPLE_DIR = @__DIR__
using DelimitedFiles
using Plots

const INPUT_PATH = joinpath(
    EXAMPLE_DIR,
    "output",
    "chang_linear_legacy_imperial_segments_uvlm_history.csv",
)
const OUTPUT_PATH = joinpath(
    EXAMPLE_DIR,
    "output",
    "chang_linear_legacy_imperial_segments_uvlm_time_histories.png",
)

isfile(INPUT_PATH) || error("Run run_chang_linear_aeroelastic.jl before plotting")

data, header = readdlm(INPUT_PATH, ',', Any, '\n'; header = true)
column_names = vec(String.(header))

function numeric_column(name::String)
    index = findfirst(==(name), column_names)
    isnothing(index) && error("Missing CSV column: $name")
    return Float64.(data[:, index])
end

time_s = numeric_column("time_s")
time_scale, time_label = maximum(time_s) < 0.1 ? (1e3, "Time (ms)") : (1.0, "Time (s)")
time = time_scale .* time_s

tip_displacement_mm = 1e3 .* numeric_column("tip_displacement_m")
tip_twist_deg = numeric_column("tip_twist_deg")
propeller_1_pitch_deg = numeric_column("propeller_1_pitch_deg")
propeller_1_yaw_deg = numeric_column("propeller_1_yaw_deg")
propeller_2_pitch_deg = numeric_column("propeller_2_pitch_deg")
propeller_2_yaw_deg = numeric_column("propeller_2_yaw_deg")
coupling_iterations = numeric_column("coupling_iterations")
state_residual = numeric_column("coupling_state_residual")
load_residual = numeric_column("coupling_load_residual")
equilibrium_residual = numeric_column("coupling_equilibrium_residual")

default(
    fontfamily = "Computer Modern",
    linewidth = 2.2,
    framestyle = :box,
    gridalpha = 0.22,
    legendfontsize = 9,
    guidefontsize = 11,
    tickfontsize = 9,
    titlefontsize = 12,
)

p_displacement = plot(
    time,
    tip_displacement_mm;
    color = :navy,
    marker = :circle,
    markersize = 3,
    xlabel = time_label,
    ylabel = "Displacement (mm)",
    title = "Wing-tip displacement",
    label = false,
)

p_twist = plot(
    time,
    tip_twist_deg;
    color = :darkorange,
    marker = :circle,
    markersize = 3,
    xlabel = time_label,
    ylabel = "Twist (deg)",
    title = "Wing-tip twist",
    label = false,
)

p_pitch = plot(
    time,
    propeller_1_pitch_deg;
    color = :royalblue,
    xlabel = time_label,
    ylabel = "Pitch (deg)",
    title = "Propeller pitch response",
    label = "Propeller 1",
    legend = :topleft,
)
plot!(p_pitch, time, propeller_2_pitch_deg; color = :crimson, linestyle = :dash, label = "Propeller 2")

p_yaw = plot(
    time,
    propeller_1_yaw_deg;
    color = :royalblue,
    xlabel = time_label,
    ylabel = "Yaw (deg)",
    title = "Propeller yaw response",
    label = "Propeller 1",
    legend = :topleft,
)
plot!(p_yaw, time, propeller_2_yaw_deg; color = :crimson, linestyle = :dash, label = "Propeller 2")

p_iterations = plot(
    time,
    coupling_iterations;
    color = :purple,
    marker = :diamond,
    markersize = 3,
    seriestype = :steppost,
    xlabel = time_label,
    ylabel = "Iterations",
    title = "Partitioned coupling iterations",
    label = false,
    ylims = (0, max(1, maximum(coupling_iterations) + 1)),
)

positive_floor = 1e-16
p_residuals = plot(
    time,
    max.(state_residual, positive_floor);
    color = :seagreen,
    yscale = :log10,
    xlabel = time_label,
    ylabel = "Residual",
    title = "Coupling convergence",
    label = "State",
    legend = :bottomright,
)
plot!(p_residuals, time, max.(load_residual, positive_floor); color = :firebrick, label = "Load")
plot!(p_residuals, time, max.(equilibrium_residual, positive_floor); color = :black, linestyle = :dot, label = "Equilibrium")

figure = plot(
    p_displacement,
    p_twist,
    p_pitch,
    p_yaw,
    p_iterations,
    p_residuals;
    layout = (3, 2),
    size = (1600, 1200),
    margin = 6Plots.mm,
    plot_title = "Chang response with legacy Imperial-segment UVLM loads",
    plot_titlefontsize = 16,
)

mkpath(dirname(OUTPUT_PATH))
savefig(figure, OUTPUT_PATH)
println("Saved time-history figure to $OUTPUT_PATH")
