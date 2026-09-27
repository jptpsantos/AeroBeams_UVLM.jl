# Summarize the diagnostic histories; does not rerun either solver.
using DelimitedFiles, LinearAlgebra, Statistics
include(joinpath(@__DIR__, "PazyWingUVLMCoupling.jl"))
const Plots = AeroBeams.Plots

function load_history(folder)
    data, header = readdlm(joinpath(folder, "history.tsv"), '\t', Float64; header=true)
    return data
end

function summarize_history(io, data, name)
    println(io, "\n", name, ": ", size(data,1), " accepted steps, last t=", data[end,1])
    println(io, "window_s | max mismatch/c | twist range_deg | max angular speed_rad/s | max nodal force_N | peak twist frequency_Hz | spectral power above 100Hz_fraction")
    for (a,b) in ((0.,0.5), (0.5,1.), (1.,1.5), (1.5,2.), (2.,2.5), (2.5,2.8), (2.8,3.01))
        ids = findall(t -> a <= t < b, data[:,1])
        length(ids) >= 20 || continue
        block = data[ids,:]
        y = block[:,3]
        n = length(y)
        t = block[:,1] .- mean(block[:,1])
        y = y - hcat(ones(n),t) * (hcat(ones(n),t) \ y)
        window = 0.5 .- 0.5 .* cos.(2pi .* (0:n-1) ./ (n-1))
        power = abs2.(AeroBeams.rfft(y .* window))
        frequency = collect(0:length(power)-1) ./ (n * median(diff(block[:,1])))
        peak = frequency[argmax(power[2:end])+1]
        fraction = sum(power[frequency .> 100]) / sum(power[2:end])
        println(io, (a,b), " | ", maximum(block[:,4]), " | ", extrema(block[:,3]),
            " | ", maximum(block[:,15]), " | ", maximum(block[:,6]), " | ", peak, " | ", fraction)
    end
    println(io, "Surface area range: ", extrema(data[:,12]))
    println(io, "Minimum wake-panel area range: ", extrema(data[:,13]))
    println(io, "Maximum surface-motion speed: ", maximum(data[:,10]))
    println(io, "Maximum wake speed: ", maximum(data[:,11]))
    println(io, "Sum of displacement-based interface work differences (NOT a full energy balance): ", sum(data[:,19]), " J")
end

function main()
    loose_folder, strong_folder = ARGS
    loose, strong = load_history(loose_folder), load_history(strong_folder)
    report_path = joinpath(loose_folder, "comparison_summary.txt")
    open(report_path, "w") do io
        summarize_history(io, loose, "Loose")
        summarize_history(io, strong, "Strong")
    end
    print(read(report_path, String))
    Plots.gr()
    panels = []
    for (col, ylabel, title, logarithmic) in (
        (2, "Tip bending [m]", "Wingtip bending", false),
        (3, "Tip twist [deg]", "Wingtip twist", false),
        (4, "Mismatch / chord", "Aerodynamic prediction versus solved geometry", true),
        (15, "Angular speed [rad/s]", "Maximum beam angular speed", true),
        (9, "max |dGamma/dt|", "Circulation-rate growth", true))
        p = Plots.plot(loose[:,1], logarithmic ? max.(loose[:,col],1e-14) : loose[:,col];
            label="Loose", color=:red, xlabel="Time [s]", ylabel, title,
            yscale=logarithmic ? :log10 : :identity, linewidth=1.2)
        Plots.plot!(p, strong[:,1], logarithmic ? max.(strong[:,col],1e-14) : strong[:,col];
            label="Strong", color=:blue, linewidth=1.2)
        push!(panels, p)
    end
    p = Plots.plot(loose[:,1], cumsum(loose[:,19]); label="Loose", color=:red,
        xlabel="Time [s]", ylabel="Cumulative difference [J]",
        title="Surface-increment work difference (diagnostic)", linewidth=1.2)
    Plots.plot!(p, strong[:,1], cumsum(strong[:,19]); label="Strong", color=:blue, linewidth=1.2)
    push!(panels,p)
    figure = Plots.plot(panels...; layout=(3,2), size=(1300,1000), margin=5Plots.mm,
        titlefontsize=10, guidefontsize=9, legendfontsize=8)
    path = joinpath(loose_folder, "coupling_comparison.png")
    Plots.savefig(figure, path)
    println("Saved ", path)
end
main()
