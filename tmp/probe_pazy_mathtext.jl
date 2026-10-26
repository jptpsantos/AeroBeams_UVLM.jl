ENV["MPLBACKEND"] = "Agg"
using Plots
with(:pyplot) do
    p = plot([30.,60.], [15.,35.];
        xlabel=raw"Airspeed, $U_\infty$ [m/s]",
        ylabel=raw"Tip OOP displacement, $u_3$ [% semispan]",
        fontfamily="Times New Roman", legend=false)
    annotate!(p, 40., 32., text(raw"(a) Tip OOP displacement, $u_3$.", 12, "Times New Roman"))
    savefig(p, joinpath(@__DIR__, "pazy_mathtext_probe.png"))
end
println("Math labels exported; restored backend: ", backend_name())
