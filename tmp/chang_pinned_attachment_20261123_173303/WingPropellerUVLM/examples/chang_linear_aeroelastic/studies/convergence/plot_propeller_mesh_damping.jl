# Run this file from VS Code or with:
# julia --project=lib/WingPropellerUVLM/examples path/to/plot_propeller_mesh_damping.jl
# Reads Mesh_case1 plus the latest root run (or Mesh_case2); no simulations run.
# Add --differences to generate separate percentage-variation figures.
if !isinteractive() && "--differences" in ARGS
    # Use the native PNG workstation for reliable serif-label measurements.
    get!(ENV,"GKSwstype","140")
end
plottingProject = normpath(joinpath(@__DIR__,"..","..",".."))
if Base.active_project() != joinpath(plottingProject,"Project.toml")
    import Pkg
    Pkg.activate(plottingProject)
end

module ChangPropellerMeshDampingPlots
using Plots,ColorSchemes,Dierckx,DelimitedFiles
include(joinpath(@__DIR__,"plot_chang_convergence_aeroelastic.jl"))
const Saved = ChangConvergenceAeroelasticPlots

"""Plot pitch and yaw damping against each propeller mesh direction.

The interpolating cubic splines (k=3, s=0) use the saved damping rates in 1/s,
not relative differences. They are evaluated only within the measured range.
"""
function plot_propeller_mesh_damping(directory=Saved.DEFAULT_DIRECTORY;
    case_names=Saved.mesh_case_names(directory),display_plots=true)
    results = Saved.load_mesh_results(directory; case_names,families=Saved.FAMILIES)
    out = joinpath(directory,"plots","propeller_mesh_damping")
    mkpath(out)
    yawColor = get(ColorSchemes.colorschemes[:darkrainbow],1.0)
    figures = Any[]
    paths = String[]
    sampled = NamedTuple[]
    for (family,name,xlabel) in (("prop_radial","spanwise","Spanwise panels per blade"),
                               ("prop_chord","chordwise","Chordwise panels per blade"))
        group = sort(filter(r -> r.family == family,results.rows); by=r -> r.x)
        panels = [r.x for r in group]
        first(panels) == 2.0 || error("The $family sweep must include the 2-panel case")
        length(panels) >= 4 && all(diff(panels) .> 0) || error("Cubic splines need four distinct panel counts")
        finePanels = collect(range(first(panels),last(panels); length=200))
        figure = plot(; xlabel,ylabel="Damping rate (1/s)",xlims=(2,last(panels)),
            xticks=sort(unique(vcat(2,collect(10:10:last(panels))))),
            tickfontsize=10,legendfontsize=10,labelfontsize=12,fontfamily="Times",
            margin=4Plots.mm,left_margin=24Plots.mm,right_margin=12Plots.mm,bottom_margin=8Plots.mm,
            grid=false,legend=:topright,framestyle=:axes,tick_direction=:in,
            size=(800,420),background_color=:white,dpi=200)
        smoothValues = Vector{Float64}[]
        for (i,label) in enumerate(("Pitch","Yaw"))
            values = [r.values[i] for r in group]
            all(isfinite,values) || error("Nonfinite $label damping in $family")
            spline = Spline1D(panels,values; k=3,s=0.0)
            all(isapprox.(spline(panels),values; rtol=1e-10,atol=1e-12)) ||
                error("Spline does not reproduce saved $label values in $family")
            smooth = spline(finePanels)
            push!(smoothValues,smooth)
            plot!(figure,finePanels,smooth; label,linewidth=3,
                lc=i == 1 ? :black : yawColor,linestyle=i == 1 ? :solid : :dash)
        end
        for j in eachindex(finePanels)
            push!(sampled,(; family,panels=finePanels[j],pitch=smoothValues[1][j],yaw=smoothValues[2][j]))
        end
        for extension in (".svg",".png")
            path = joinpath(out,"propeller_$(name)_damping_spline" * extension)
            savefig(figure,path)
            push!(paths,path)
        end
        push!(figures,figure)
        display_plots && display(figure)
    end
    open(joinpath(out,"damping_mesh_samples.csv"),"w") do io
        writedlm(io,permutedims(["family","panels","pitch_lambda_per_s","yaw_lambda_per_s",
            "valid","status","reason","sources"]),',')
        for r in results.rows
            writedlm(io,permutedims(Any[r.family,r.x,r.values[1],r.values[2],
                r.valid,r.status,r.reason,r.sources]),',')
        end
    end
    open(joinpath(out,"damping_mesh_interpolated.csv"),"w") do io
        writedlm(io,permutedims(["family","panels","pitch_lambda_per_s","yaw_lambda_per_s"]),',')
        for r in sampled
            writedlm(io,permutedims(Any[r.family,r.panels,r.pitch,r.yaw]),',')
        end
    end
    open(joinpath(out,"plot_notes.md"),"w") do io
        println(io,"# Propeller mesh damping\n")
        println(io,"Sources: $(join(case_names," + ")); . denotes the latest root run. Full study settings, fingerprints and duplicate configurations are checked before merging.")
        println(io,"Two figures: spanwise (radial) panels and chordwise panels per blade. Both show pitch and yaw damping rate lambda in 1/s; negative values mean decay.")
        println(io,"Dierckx interpolating cubic splines use k=3 and s=0, evaluated on 200 points from 2 to 40 panels. All original knots are reproduced to numerical precision; no extrapolation is used.")
        println(io,"The original mesh levels are 2, 5, 7, 10, 20, 30 and 40. Black solid: pitch. Darkrainbow endpoint dashed: yaw. Linewidth 3, tick/legend font 10, axis label font 12, no grid.")
        println(io,"Smoothing is for presentation; it adds no solved mesh levels and does not change fit validity. Original values, status and sources are retained in damping_mesh_samples.csv.")
        println(io,"$(count(r -> !r.valid,results.rows)) of $(length(results.rows)) mesh points have unvalidated damping fits; these remain diagnostic data.")
    end
    foreach(println,paths)
    return (; figures,paths,results,sampled)
end

"""Cubic Hermite interpolation that preserves the shape between mesh samples.

Weighted harmonic interior slopes prevent added extrema; endpoints use secants.
"""
function interpolate_mesh_variation(panels,values,requestedPanels)
    n = length(panels)
    n == length(values) && n >= 2 || error("Need matching mesh and value arrays")
    widths = diff(panels)
    all(>(0),widths) && all(isfinite,values) || error("Invalid interpolation samples")
    all(x -> first(panels) <= x <= last(panels),requestedPanels) || error("Extrapolation is disabled")
    secants = diff(values) ./ widths
    slopes = zeros(n)
    slopes[1],slopes[end] = first(secants),last(secants)
    for j in 2:n-1
        before,after = secants[j-1],secants[j]
        if !iszero(before) && !iszero(after) && sign(before) == sign(after)
            w1 = 2 * widths[j] + widths[j-1]
            w2 = widths[j] + 2 * widths[j-1]
            slopes[j] = (w1 + w2) / (w1 / before + w2 / after)
        end
    end
    return map(requestedPanels) do x
        j = min(searchsortedlast(panels,x),n-1)
        h = widths[j]
        t = (x - panels[j]) / h
        (2 * t^3 - 3 * t^2 + 1) * values[j] +
        (t^3 - 2 * t^2 + t) * h * slopes[j] +
        (-2 * t^3 + 3 * t^2) * values[j+1] +
        (t^3 - t^2) * h * slopes[j+1]
    end
end

"""Plot absolute percentage changes between successive solved mesh levels.

Differences use the current (finer) damping magnitude as denominator. Smoothing
is applied after computing the differences, never to generate refinement data.
"""
function plot_propeller_mesh_differences(directory=Saved.DEFAULT_DIRECTORY;
    case_names=Saved.mesh_case_names(directory),display_plots=true)
    results = Saved.load_mesh_results(directory; case_names,families=Saved.FAMILIES)
    out = joinpath(directory,"plots","propeller_mesh_damping")
    mkpath(out)
    yawColor = get(ColorSchemes.colorschemes[:darkrainbow],1.0)
    figures = Any[]
    paths = String[]
    sampled = NamedTuple[]
    changes = NamedTuple[]
    for (family,name,xlabel) in (("prop_radial","spanwise","Spanwise panels per blade"),
                               ("prop_chord","chordwise","Chordwise panels per blade"))
        group = sort(filter(r -> r.family == family,results.rows); by=r -> r.x)
        points = group[2:end]
        panels = [r.x for r in points]
        length(panels) >= 4 || error("Cubic splines need four refinement comparisons")
        finePanels = collect(range(first(panels),last(panels); length=200))
        figure = plot(; xlabel,ylabel="Damping variation (%)",xlims=(2,last(panels)),
            xticks=sort(unique(vcat(2,first(panels),collect(10:10:last(panels))))),
            tickfontsize=10,legendfontsize=10,labelfontsize=12,fontfamily="Computer Modern",
            margin=4Plots.mm,left_margin=32Plots.mm,right_margin=12Plots.mm,bottom_margin=8Plots.mm,
            grid=false,legend=:topright,framestyle=:axes,tick_direction=:in,
            size=(800,420),background_color=:white,dpi=200)
        smoothValues = Vector{Float64}[]
        for (i,label) in enumerate(("Pitch","Yaw"))
            values = [r.successive_percent[i] for r in points]
            all(x -> isfinite(x) && x >= 0,values) || error("Undefined $label variation in $family")
            knots = interpolate_mesh_variation(panels,values,panels)
            all(isapprox.(knots,values; rtol=1e-10,atol=1e-12)) || error("Variation spline misses saved knots")
            smooth = interpolate_mesh_variation(panels,values,finePanels)
            all(x -> isfinite(x) && x >= 0,smooth) || error("Invalid interpolated variation")
            push!(smoothValues,smooth)
            color = i == 1 ? :black : yawColor
            plot!(figure,finePanels,smooth; label,linewidth=3,lc=color,
                linestyle=i == 1 ? :solid : :dash)
            scatter!(figure,panels,values; label="",color,markersize=3,markerstrokewidth=0)
        end
        ymax = 1.08 * maximum(maximum, smoothValues)
        plot!(figure; ylims=(0.0,max(ymax,1e-3)))
        for j in eachindex(finePanels)
            push!(sampled,(; family,panels=finePanels[j],pitch=smoothValues[1][j],yaw=smoothValues[2][j]))
        end
        for (j,r) in enumerate(group)
            previousPanels = j == 1 ? NaN : group[j-1].x
            push!(changes,(; family,previous_panels=previousPanels,panels=r.x,
                pitch=r.successive_percent[1],yaw=r.successive_percent[2],
                valid=r.successive_valid,sources=r.sources))
        end
        for extension in (".svg",".png")
            path = joinpath(out,"propeller_$(name)_damping_variation_percent" * extension)
            savefig(figure,path)
            push!(paths,path)
        end
        push!(figures,figure)
        display_plots && display(figure)
    end
    open(joinpath(out,"damping_variation_samples.csv"),"w") do io
        writedlm(io,permutedims(["family","previous_panels","panels","pitch_variation_percent",
            "yaw_variation_percent","comparison_valid","current_sources"]),',')
        for r in changes
            writedlm(io,permutedims(Any[r.family,r.previous_panels,r.panels,r.pitch,r.yaw,r.valid,r.sources]),',')
        end
    end
    open(joinpath(out,"damping_variation_interpolated.csv"),"w") do io
        writedlm(io,permutedims(["family","panels","pitch_variation_percent","yaw_variation_percent"]),',')
        for r in sampled
            writedlm(io,permutedims(Any[r.family,r.panels,r.pitch,r.yaw]),',')
        end
    end
    open(joinpath(out,"damping_variation_notes.md"),"w") do io
        println(io,"# Percentage variation between propeller mesh refinements\n")
        println(io,"Sources: $(join(case_names," + ")); . denotes the latest root run. Existing damping figures and their data files are preserved.\n")
        println(io,"Variation (%) = 100 * abs(lambda_i - lambda_previous) / abs(lambda_i). The denominator is the magnitude of the current, finer-mesh damping rate, consistent with the earlier percentage plots.")
        println(io,"Comparisons are 2 -> 5, 5 -> 7, 7 -> 10, 10 -> 20, 20 -> 30 and 30 -> 40 panels. Each value is placed at the finer panel count. The 2-panel point has no previous result: NaN in the CSV and omitted from the curves.\n")
        println(io,"Percentages are computed from solved mesh values before interpolation. Dots mark the actual comparisons. Shape-preserving cubic Hermite interpolation uses weighted harmonic interior slopes and secant endpoint slopes on 200 points. Each segment stays between its measured endpoint values, avoiding artificial extrema or negative variations. No extrapolation is applied. The original damping-rate figures retain their Dierckx cubic splines.")
        println(io,"Interpolation is only for presentation. Fit validity is retained in damping_variation_samples.csv; small variation alone does not validate a rejected damping fit.")
        println(io,"The ordinate measures changes between adjacent tested levels, not the error against an exact solution; refinement increments are unequal.")
    end
    foreach(println,paths)
    return (; figures,paths,results,sampled,changes)
end
end

if abspath(PROGRAM_FILE) == (@__FILE__) || (isinteractive() && (@__MODULE__) === Main)
    differences = "--differences" in ARGS
    directories = filter(!=("--differences"),ARGS)
    directory = isempty(directories) ? ChangPropellerMeshDampingPlots.Saved.DEFAULT_DIRECTORY : only(directories)
    if differences
        propeller_mesh_differences = ChangPropellerMeshDampingPlots.plot_propeller_mesh_differences(
            directory; display_plots=isinteractive())
    else
        propeller_mesh_damping = ChangPropellerMeshDampingPlots.plot_propeller_mesh_damping(directory)
    end
end
