# Plot the selected independent upsweep and pulsed downsweep together.
# This reads saved histories only; no simulation is executed.
include(joinpath(@__DIR__, "plot_pazy_LCO.jl"))
using SHA

function plot_pazy_LCO_comparison(;
    upsweep_csv=joinpath(@__DIR__, "output", "LCO_independent_sweeps", "pazy_linear_LCO_upsweep.csv"),
    downsweep_csv=joinpath(@__DIR__, "output", "LCO_downsweep_pulse", "pazy_linear_LCO_downsweep.csv"),
    output_directory=joinpath(@__DIR__, "output", "LCO_selected_comparison"))
    branches = []
    sources = String[]
    for (path, label, color) in ((upsweep_csv, "Upsweep", :blue),
                                (downsweep_csv, "Downsweep", :red))
        bytes = read(path)
        data, header = readdlm(IOBuffer(bytes), ',', Float64; header=true)
        column(name) = data[:,something(findfirst(==(name), vec(header)))]
        speed = column("airspeed_mps")
        strain = column("root_longitudinal_strain_ue")
        time = column("time_s")
        all(isfinite, speed) && all(isfinite, strain) || error("Nonfinite data in $path")
        all(diff(time) .> 0) || error("Non-increasing time in $path")
        all((label == "Upsweep" ? (>) : (<)).(diff(speed), 0)) ||
            error("Unexpected speed direction in $path")
        push!(branches, (; airspeed=speed, strain, label, color))
        push!(sources, "$label: $(abspath(path))\nSHA256: $(bytes2hex(sha256(bytes)))\n" *
            "Samples: $(length(time)); time: $(first(time)) -> $(last(time)) s; " *
            "speed: $(first(speed)) -> $(last(speed)) m/s")
    end
    figure = pazy_LCO_strain_plot(branches)
    mkpath(output_directory)
    output_path = joinpath(output_directory, "pazy_LCO_selected_strain_vs_airspeed.png")
    Plots.savefig(figure, output_path)
    write(joinpath(output_directory, "sources.txt"),
        "Instantaneous saved strain versus speed; original sample order; no filtering or rescaling.\n\n" *
        join(sources, "\n\n") * "\n")
    println(join(sources, "\n\n"))
    println("Saved plot: ", output_path)
    return figure
end

if abspath(PROGRAM_FILE) == @__FILE__
    plot_pazy_LCO_comparison()
end
