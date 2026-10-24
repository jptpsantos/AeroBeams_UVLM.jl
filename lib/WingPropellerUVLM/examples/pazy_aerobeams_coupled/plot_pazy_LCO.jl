# Plot saved LCO strain histories without running either simulation.
using DelimitedFiles
import Plots

function pazy_lco_branch_names(selection; option_name="plot_branches")
    selection == :both && return ("upsweep", "downsweep")
    selection == :upsweep && return ("upsweep",)
    selection == :downsweep && return ("downsweep",)
    throw(ArgumentError("$option_name must be :both, :upsweep, or :downsweep"))
end

"""Plot instantaneous strain against speed, preserving each branch's time order."""
function pazy_LCO_strain_plot(branches; xlims=:auto, ylims=:auto)
    figure = Plots.plot(;
        xlabel="Airspeed [m/s]", ylabel="Root longitudinal strain [με]",
        fontfamily="Computer Modern", guidefontsize=16, tickfontsize=12,
        legendfontsize=12, legend=:bottomleft,
        background_color_legend=:white, foreground_color_legend=:gray,
        framestyle=:axes, grid=true, gridalpha=0.2,
        size=(800, 550), dpi=300, margin=5Plots.mm, xlims, ylims)
    for branch in branches
        Plots.plot!(figure, branch.airspeed, branch.strain;
            label=branch.label, color=branch.color, linewidth=0.5, linealpha=0.65)
    end
    return figure
end

"""Regenerate the strain-versus-speed PNG from available branch CSVs.

Both branches are used when available. A missing branch is reported and omitted,
so a completed upsweep can be inspected while the downsweep is still running.
The default reads the independent runs. Use `sweep_mode=:continuous` to inspect
the continuous test, or pass an explicit `output_directory`.
Optional `xlims` and `ylims` set axis limits; automatic limits retain all data.
`plot_branches` selects :both (default), :upsweep, or :downsweep. Single-branch
figures have a branch suffix so they do not replace the combined figure.
"""
function plot_pazy_LCO(;
    sweep_mode=:independent,
    output_directory=joinpath(@__DIR__, "output",
        sweep_mode == :continuous ? "LCO_continuous_sweeps" : "LCO_independent_sweeps"),
    plot_branches=:both, xlims=:auto, ylims=:auto)
    @assert sweep_mode in (:continuous, :independent)
    selected_names = pazy_lco_branch_names(plot_branches)
    branches = []
    for (name, label, color) in (("upsweep", "Upsweep", :blue),
                                ("downsweep", "Downsweep", :red))
        name in selected_names || continue
        path = joinpath(output_directory, "pazy_linear_LCO_$(name).csv")
        if !isfile(path)
            @warn "Branch CSV not available; omitting $label" path
            continue
        end
        data, header = readdlm(path, ',', Float64; header=true)
        speed_column = findfirst(==("airspeed_mps"), vec(header))
        strain_column = findfirst(==("root_longitudinal_strain_ue"), vec(header))
        if isnothing(speed_column) || isnothing(strain_column)
            throw(ArgumentError("Missing airspeed or root strain column in $path"))
        end
        push!(branches, (; airspeed=data[:, speed_column],
            strain=data[:, strain_column], label, color))
    end
    isempty(branches) && throw(ArgumentError("No LCO branch CSVs found in $output_directory"))
    figure = pazy_LCO_strain_plot(branches; xlims, ylims)
    suffix = plot_branches == :both ? "" : "_$(plot_branches)"
    output_path = joinpath(output_directory, "pazy_linear_LCO_root_strain_vs_airspeed$suffix.png")
    Plots.savefig(figure, output_path)
    println("Strain-versus-airspeed plot saved in: $output_path")
    return figure
end

#if abspath(PROGRAM_FILE) == @__FILE__
    plot_pazy_LCO()
#end
