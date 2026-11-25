# Read one saved aeroelastic campaign and plot its sampled growth rates and
# successive differences in the VFS plot_sensitivity style.
using DelimitedFiles,Printf,TOML
get!(ENV,"GKSwstype","100")
using Plots
include(joinpath(@__DIR__,"plot_chang_convergence_aeroelastic.jl"))
const Saved = ChangConvergenceAeroelasticPlots
const WING_METRICS = ("wing_bending_lambda_per_s","wing_torsion_lambda_per_s",
    "wing_bending_frequency_hz","wing_torsion_frequency_hz")
const WING_SIGNALS = ("Wing bending","Wing torsion")
const PROPELLER_SIGNALS = ("Propeller pitch","Propeller yaw")

const AXIS_LABELS = (
    "wing_span" => "Wing spanwise panels",
    "wing_chord" => "Wing chordwise panels",
    "prop_radial" => "Radial panels per blade",
    "prop_chord" => "Chordwise panels per blade",
    "wake_revolutions" => "Wake cutoff (rev.)",
    "azimuth_deg" => "Azimuthal step (deg.)"
)

# Same shape-preserving cubic interpolation as plot_sensitivity_percent.jl.
function interpolate_sensitivity(x,y,query)
    n = length(x)
    n == length(y) && n >= 2 || error("Need at least two matched samples")
    widths = diff(x)
    all(>(0),widths) && all(isfinite,y) || error("Invalid interpolation samples")
    all(q -> first(x) <= q <= last(x),query) || error("Interpolation cannot extrapolate")
    secants = diff(y) ./ widths
    slopes = zeros(n)
    if n == 2
        slopes .= first(secants)
    else
        for i in 2:n-1
            if secants[i-1] * secants[i] > 0
                w1,w2 = 2widths[i] + widths[i-1],widths[i] + 2widths[i-1]
                slopes[i] = (w1 + w2) / (w1 / secants[i-1] + w2 / secants[i])
            end
        end
        function endpoint(h1,h2,d1,d2)
            d = ((2h1 + h2) * d1 - h1 * d2) / (h1 + h2)
            sign(d) == sign(d1) || return 0.0
            return sign(d1) != sign(d2) && abs(d) > 3abs(d1) ? 3d1 : d
        end
        slopes[1] = endpoint(widths[1],widths[2],secants[1],secants[2])
        slopes[end] = endpoint(widths[end],widths[end-1],secants[end],secants[end-1])
    end
    return map(query) do q
        i = min(searchsortedlast(x,q),n-1)
        h = widths[i]
        t = (q - x[i]) / h
        (2t^3 - 3t^2 + 1) * y[i] + (t^3 - 2t^2 + t) * h * slopes[i] +
        (-2t^3 + 3t^2) * y[i+1] + (t^3 - t^2) * h * slopes[i+1]
    end
end

function sensitivity_tick_step(span)
    rawStep = span / 6
    magnitude = 10.0^floor(log10(rawStep))
    candidates = (1.0,2.0,2.5,5.0,10.0)
    return magnitude * argmin(step -> abs(step - rawStep/magnitude),candidates)
end

function sensitivity_panel(group,family,xlabel,scale,signals;showLegend)
    growth = scale == :growth
    sampledRows = growth ? group : group[2:end]
    rows = sort(sampledRows;by=r -> r.x)
    x = [r.x for r in rows]
    query = range(first(x),last(x);length=600)
    field = growth ? :values : scale == :percent ? :successive_percent : :previous
    smooth = map(1:2) do i
        y = [getproperty(r,field)[i] for r in rows]
        knots = interpolate_sensitivity(x,y,x)
        all(isapprox.(knots,y;atol=1e-12,rtol=1e-10)) || error("Interpolation missed a solved point")
        result = interpolate_sensitivity(x,y,query)
        minimum(y)-1e-12 <= minimum(result) <= maximum(result) <= maximum(y)+1e-12 ||
            error("Interpolation overshot the solved data")
        return result
    end
    if growth
        lo,hi = extrema(Iterators.flatten(smooth))
        pad = 0.08max(hi-lo,1e-3)
        step = sensitivity_tick_step(hi-lo+2pad)
        lower = step * floor((lo-pad)/step)
        upper = step * ceil((hi+pad)/step)
    else
        headroom = showLegend ? 1.32 : 1.08
        peak = headroom * max(maximum(maximum,smooth),1e-3)
        step = sensitivity_tick_step(peak)
        lower,upper = 0.0,step * ceil(peak/step)
    end
    xstep = family == "azimuth_deg" ? 5 : family == "wake_revolutions" ? 1 : 10
    xmax = xstep * ceil(maximum(r.x for r in group) / xstep)
    azimuthDifferences = family == "azimuth_deg" && !growth
    if azimuthDifferences
        padding = 0.04 * (last(x) - first(x))
        limits = (first(x)-padding,last(x)+padding)
        ticks = (x,[@sprintf("%g",stepValue) for stepValue in x])
        xlabel = "Finer azimuthal step (deg.)"
    else
        limits = (0,xmax)
        ticks = 0:xstep:xmax
    end
    ylabel = growth ? "Growth rate λ (s⁻¹)" : scale == :percent ?
        "Difference between refinements (%)" : "Damping difference (s⁻¹)"
    p = plot(query,smooth[1];label=signals[1],color=:black,linestyle=:solid,linewidth=3,
        fontfamily="sans-serif",tickfontsize=10,legendfontsize=10,labelfontsize=12,
        margin=4Plots.mm,grid=false,legend=showLegend ? (azimuthDifferences ? :topleft : :topright) : false,
        foreground_color_legend=nothing,background_color_legend=nothing,
        framestyle=:axes,mirror=false,tick_direction=:in,
        xlabel,ylabel,xlims=limits,xticks=ticks,xflip=azimuthDifferences,
        yscale=:identity,ylims=(lower,upper),yticks=lower:step:upper,size=(600,400),dpi=300)
    plot!(p,query,smooth[2];label=signals[2],color=:red,linestyle=:dash,linewidth=3)
    if azimuthDifferences
        scatter!(p,x,[getproperty(r,field)[1] for r in rows];
            label="",markercolor=:white,markerstrokecolor=:black,markersize=4)
        scatter!(p,x,[getproperty(r,field)[2] for r in rows];
            label="",markercolor=:red,markerstrokewidth=0,markersize=2.2)
    end
    return p
end

function plot_saved_sensitivity(directory)
    target = get(TOML.parsefile(joinpath(directory,"study.toml"))["settings"]["damping"],"target","propeller")
    metrics = target == "wing" ? WING_METRICS : Saved.METRICS
    signals = target == "wing" ? WING_SIGNALS : PROPELLER_SIGNALS
    data = Saved.load_current_results(directory;metrics)
    groups = [(;family,xlabel,signals,rows=filter(r -> r.family==family,data.rows)) for (family,xlabel) in AXIS_LABELS]
    return save_sensitivity_plots(groups,joinpath(directory,"plots","sensitivity_style"),
        data.settings["physical"]["speed_mps"])
end

function save_sensitivity_plots(groups,output,airspeed)
    filter!(g -> !isempty(g.rows),groups)
    isempty(groups) && error("No supported sweep families in the saved run")
    all(g -> length(g.rows)>=3,groups) || error("Each plotted sweep needs at least three solved levels")
    mkpath(output)
    speed = @sprintf("%g m/s",airspeed)
    flag = all(g -> all(r -> r.valid,g.rows),groups) ? "" : " (diagnostic fits present)"
    for (scale,title) in (
        (:growth,"Fitted growth rates"),
        (:absolute,"Absolute successive differences"),
        (:percent,"Percentage successive differences")
    )
        figures = Plots.Plot[]
        legends = Set{Tuple{String,String}}()
        for g in groups
            p = sensitivity_panel(g.rows,g.family,g.xlabel,scale,g.signals;showLegend=g.signals ∉ legends)
            push!(legends,g.signals)
            name = "sensitivity_$(g.family)_$scale"
            for ext in ("png","svg")
                savefig(p,joinpath(output,"$name.$ext"))
            end
            push!(figures,p)
        end
        nrows = cld(length(figures),2)
        combined = plot(figures...;layout=(nrows,2),size=(1200,400nrows),dpi=150,
            plot_title="$speed · $title$flag",plot_titlefontsize=14)
        savefig(combined,joinpath(output,"convergence_$scale.png"))
    end
    azimuth = findfirst(g -> g.family == "azimuth_deg",groups)
    if !isnothing(azimuth)
        g = groups[azimuth]
        percent = sensitivity_panel(g.rows,g.family,g.xlabel,:percent,g.signals;showLegend=true)
        absolute = sensitivity_panel(g.rows,g.family,g.xlabel,:absolute,g.signals;showLegend=false)
        plot!(absolute;left_margin=8Plots.mm,bottom_margin=8Plots.mm)
        plot!(percent;left_margin=8Plots.mm,bottom_margin=8Plots.mm)
        azimuthFigure = plot(absolute,percent;layout=(1,2),size=(1200,480),dpi=150,
            plot_title="$speed · Azimuthal step refinements$flag",plot_titlefontsize=14)
        for ext in ("png","svg")
            savefig(azimuthFigure,joinpath(output,"azimuth_step_differences.$ext"))
        end
    end
    open(joinpath(output,"successive_differences.csv"),"w") do io
        println(io,"family,previous_level,refined_level,response,lambda_previous_per_s,lambda_refined_per_s,absolute_difference_per_s,percentage_difference,valid_pair")
        for group in groups, j in 2:length(group.rows), i in 1:2
            previous,current = group.rows[j-1],group.rows[j]
            fields = (group.family,previous.x,current.x,group.signals[i],previous.values[i],current.values[i],
                current.previous[i],current.successive_percent[i],current.successive_valid)
            writedlm(io,permutedims(collect(fields)),',')
        end
    end
    return output
end

if abspath(PROGRAM_FILE) == @__FILE__
    println(plot_saved_sensitivity(abspath(only(ARGS))))
end
