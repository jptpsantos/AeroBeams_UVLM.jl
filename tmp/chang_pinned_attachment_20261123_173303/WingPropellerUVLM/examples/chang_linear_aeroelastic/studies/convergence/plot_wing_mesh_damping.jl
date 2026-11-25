# Run from VS Code or with:
# julia --project=lib/WingPropellerUVLM/examples path/to/plot_wing_mesh_damping.jl
# Add --differences for successive percentage variations, as for the propeller.
# Only saved tables are read; no simulations or damping fits are rerun.
if !isinteractive()
    get!(ENV,"GKSwstype","140")
end

module ChangWingMeshDampingPlots
include(joinpath(@__DIR__,"plot_propeller_mesh_damping.jl"))
using Plots,ColorSchemes,Dierckx,DelimitedFiles,Printf
const Propeller = ChangPropellerMeshDampingPlots
const Saved = Propeller.Saved
const DEFAULT_DIRECTORY = normpath(joinpath(@__DIR__,"..","..","output","convergence_aeroelastic_wing"))
const METRICS = ("wing_bending_lambda_per_s","wing_torsion_lambda_per_s",
    "wing_bending_frequency_hz","wing_torsion_frequency_hz")

"""Plot saved wing damping or its successive percentage variation.

Uses the propeller colors, line styles and interpolation procedures. Damping
rates use interpolating Dierckx cubic splines. Percentage changes are computed
on solved levels first, then use shape-preserving cubic Hermite interpolation.
Rejected fits retain their status and remain diagnostic data.
"""
function plot_wing_mesh_damping(directory=DEFAULT_DIRECTORY;
    differences=false,display_plots=isinteractive())
    directory = abspath(directory)
    results = Saved.load_current_results(directory; metrics=METRICS)
    out = joinpath(directory,"plots","wing_mesh_damping")
    mkpath(out)
    torsionColor = get(ColorSchemes.colorschemes[:darkrainbow],1.0)
    figures = Any[]
    paths = String[]
    sampled = NamedTuple[]
    changes = NamedTuple[]
    suffix = differences ? "variation_percent" : "spline"
    for (family,name,xlabel) in (("wing_span","spanwise","Spanwise panels"),
                               ("wing_chord","chordwise","Chordwise panels"))
        group = sort(filter(r -> r.family == family,results.rows); by=r -> r.x)
        length(group) >= 4 || error("$family needs four distinct panel counts")
        all(diff([r.x for r in group]) .> 0) || error("Repeated panel counts in $family")
        points = differences ? group[2:end] : group
        panels = [r.x for r in points]
        finePanels = collect(range(first(panels),last(panels); length=200))
        ticks = sort(unique(vcat(first(group).x,first(panels),collect(10:10:last(panels)))))
        figure = plot(; xlabel,
            ylabel=differences ? "Damping variation (%)" : "Damping rate (1/s)",
            xlims=(first(group).x,last(panels)),xticks=ticks,
            tickfontsize=10,legendfontsize=10,labelfontsize=12,
            fontfamily=differences ? "Computer Modern" : "Times",
            margin=4Plots.mm,left_margin=(differences ? 32 : 24) * Plots.mm,
            right_margin=12Plots.mm,bottom_margin=8Plots.mm,
            grid=false,legend=!differences && name == "chordwise" ? :topleft : :topright,
            framestyle=:axes,tick_direction=:in,
            size=(800,420),background_color=:white,dpi=200)
        smoothValues = Vector{Float64}[]
        for (i,label) in enumerate(("Bending","Torsion"))
            values = [differences ? r.successive_percent[i] : r.values[i] for r in points]
            all(isfinite,values) || error("Undefined $label damping or variation in $family")
            interpolant = differences ?
                x -> Propeller.interpolate_mesh_variation(panels,values,x) :
                Spline1D(panels,values; k=3,s=0.0)
            all(isapprox.(interpolant(panels),values; rtol=1e-10,atol=1e-12)) ||
                error("Interpolation misses saved $label values in $family")
            smooth = interpolant(finePanels)
            all(isfinite,smooth) || error("Nonfinite interpolated damping")
            differences && !all(>=(0),smooth) && error("Negative interpolated variation")
            push!(smoothValues,smooth)
            color = i == 1 ? :black : torsionColor
            plot!(figure,finePanels,smooth; label,linewidth=3,lc=color,
                linestyle=i == 1 ? :solid : :dash)
            if differences
                scatter!(figure,panels,values; label="",color,markersize=3,markerstrokewidth=0)
            end
        end
        if differences
            ymax = 1.08 * maximum(maximum,smoothValues)
            plot!(figure; ylims=(0.0,max(ymax,1e-3)))
        end
        for j in eachindex(finePanels)
            push!(sampled,(; family,panels=finePanels[j],
                bending=smoothValues[1][j],torsion=smoothValues[2][j]))
        end
        for (j,r) in enumerate(group)
            push!(changes,(; family,previous_panels=j == 1 ? NaN : group[j-1].x,
                panels=r.x,bending=r.successive_percent[1],torsion=r.successive_percent[2],
                valid=r.successive_valid,sources=r.sources))
        end
        for extension in (".svg",".png")
            path = joinpath(out,"wing_$(name)_damping_$(suffix)" * extension)
            savefig(figure,path)
            push!(paths,path)
        end
        push!(figures,figure)
        display_plots && display(figure)
    end
    dataName = differences ? "damping_variation" : "damping_mesh"
    metricNames = differences ? ("bending_variation_percent","torsion_variation_percent") :
        ("wing_bending_lambda_per_s","wing_torsion_lambda_per_s")
    open(joinpath(out,dataName * "_samples.csv"),"w") do io
        if differences
            writedlm(io,permutedims(["family","previous_panels","panels",metricNames...,
                "comparison_valid","current_sources"]),',')
            for r in changes
                writedlm(io,permutedims(Any[r.family,r.previous_panels,r.panels,
                    r.bending,r.torsion,r.valid,r.sources]),',')
            end
        else
            writedlm(io,permutedims(["family","panels",metricNames...,
                "valid","status","reason","sources"]),',')
            for r in results.rows
                writedlm(io,permutedims(Any[r.family,r.x,r.values[1],r.values[2],
                    r.valid,r.status,r.reason,r.sources]),',')
            end
        end
    end
    open(joinpath(out,dataName * "_interpolated.csv"),"w") do io
        writedlm(io,permutedims(["family","panels",metricNames...]),',')
        for r in sampled
            writedlm(io,permutedims(Any[r.family,r.panels,r.bending,r.torsion]),',')
        end
    end
    open(joinpath(out,differences ? "damping_variation_notes.md" : "plot_notes.md"),"w") do io
        println(io,"# Wing mesh damping",differences ? " percentage variation\n" : "\n")
        println(io,"Source: $directory. Reads convergence.csv, case_matrix.csv and study.toml; no fits or simulations are rerun.\n")
        println(io,"Black solid: bending. Darkrainbow endpoint dashed: torsion. Linewidth 3, tick/legend font 10, axis labels 12, no grid; styling follows plot_propeller_mesh_damping.jl.\n")
        if differences
            println(io,"Variation (%) = 100 * abs(lambda_i - lambda_previous) / abs(lambda_i). Each point is placed at the finer panel count. The first level has no predecessor: NaN in the CSV, omitted from the curves.")
            println(io,"Percentages are computed before interpolation. Dots mark solved comparisons. Shape-preserving cubic Hermite interpolation uses the propeller helper, weighted harmonic interior slopes and secant endpoint slopes, evaluated on 200 points with no extrapolation.")
            println(io,"These are adjacent-level changes, not errors against an exact solution; refinement increments are unequal. Small variation does not validate a rejected fit.\n")
        else
            println(io,"Damping rates are lambda in 1/s; negative means decay. Dierckx cubic splines (k=3, s=0) reproduce all saved knots and use 200 points within each measured range, with no extrapolation. Interpolation is for presentation and adds no solved mesh levels.\n")
        end
        for family in ("wing_span","wing_chord")
            group = filter(r -> r.family == family,results.rows)
            println(io,"$family panel counts: ",join([@sprintf("%g",r.x) for r in group],", "),".")
        end
        println(io,"\n$(count(r -> !r.valid,results.rows)) of $(length(results.rows)) saved fits are unvalidated. These figures remain diagnostic; original validity/status/reasons are retained in damping_mesh_samples.csv, and comparison validity in damping_variation_samples.csv.")
        println(io,"Span refinement changes both structural and aerodynamic meshes. The chord sweep keeps the span count fixed. This study does not establish a flutter boundary.")
    end
    foreach(println,paths)
    return (; figures,paths,results,sampled,changes)
end
end

if abspath(PROGRAM_FILE) == (@__FILE__) || isinteractive()
    differences = "--differences" in ARGS
    directories = filter(!=("--differences"),ARGS)
    directory = isempty(directories) ? ChangWingMeshDampingPlots.DEFAULT_DIRECTORY : only(directories)
    wing_mesh_damping = ChangWingMeshDampingPlots.plot_wing_mesh_damping(directory; differences)
end
