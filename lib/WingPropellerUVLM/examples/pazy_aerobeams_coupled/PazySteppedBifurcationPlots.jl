# CSV-only plotting; including this file never loads an aeroelastic solver.
module PazySteppedBifurcationPlots
using Plots, DelimitedFiles
include(joinpath(@__DIR__,"plot_pazy_bifurcation_interpolated.jl"))
import .PazyInterpolatedBifurcation as Interpolation
const Paper=Interpolation.Paper
export plot_stepped_bifurcation

"""Plot ascending/descending envelopes independently, preserving jump direction.
Inputs are the branch CSVs under a stepped run directory. Every invocation writes
to a fresh figure directory; raw CSVs and previous figures are never overwritten.
Solid/open = upsweep; dashed/filled = downsweep. Unverified knots have black crosses.
"""
function plot_stepped_bifurcation(run_directory; linewidth=2.8,markersize=3.4,
    formats=("png","pdf"),save_panels=true,points_per_interval=32,
    figure_size=(1300,1040),panel_size=(650,460),dpi=300,xlims=nothing,
    ylims=(nothing,nothing,nothing,nothing),show_local_extrema=true)
    all(v -> isfinite(v) && v > 0,(linewidth,markersize,dpi)) || throw(ArgumentError("Invalid positive plot dimensions"))
    length(ylims) == 4 || throw(ArgumentError("ylims requires four entries"))
    !isempty(formats) && all(f -> f in ("png","pdf","svg"),formats) || throw(ArgumentError("Invalid figure formats"))
    cases=NamedTuple[]
    for direction in ("upsweep","downsweep")
        path=joinpath(abspath(run_directory),direction,"bifurcation_summary.csv")
        isfile(path) || continue
        data,_=Paper.read_table(path)
        size(data,1) == 0 && continue
        angle=Paper.saved_incidence(path)
        dataset=Interpolation.load_case(path,angle;points_per_interval,max_gap=nothing)
        push!(cases,(;direction,dataset))
    end
    isempty(cases) && throw(ArgumentError("No completed branch summaries under $run_directory"))
    root=joinpath(abspath(run_directory),"figures")
    mkpath(root)
    directory=mktempdir(root;prefix="plots_",cleanup=false)
    previous_backend=Plots.backend_name()
    try
        Plots.gr()
        panels,captions,saved=Any[],Any[],String[]
        common=(;fontfamily="Times",guidefontsize=13,tickfontsize=11,legendfontsize=10,
            framestyle=:box,grid=true,gridalpha=0.18,foreground_color_grid=:gray,
            tick_direction=:in,background_color=:white,background_color_legend=:white,
            left_margin=6Plots.mm,right_margin=4Plots.mm,bottom_margin=4Plots.mm,top_margin=3Plots.mm)
        for (i,response) in enumerate(Paper.PANEL_RESPONSES)
            p=plot(;common...,xlabel=Paper.NativeTeXLabel("Airspeed, U_{∞} [m/s]"),
                ylabel=Paper.PANEL_LABELS[i],legend=i == 1 ? :topleft : false)
            for case in cases
                d=case.dataset.panels[i]
                color=case.dataset.color
                up=case.direction == "upsweep"
                fill=up ? :white : color
                if show_local_extrema
                    for (kind,marker) in (("minimum",:dtriangle),("maximum",:utriangle))
                        selected=d.peak_kind .== kind
                        scatter!(p,d.peak_speed[selected],d.peak_value[selected];marker,markersize,
                            markercolor=fill,markerstrokecolor=color,markerstrokewidth=0.6,markeralpha=0.5,label=false)
                    end
                end
                for curve in d.curves
                    for values in (curve.minimum,curve.maximum)
                        plot!(p,curve.speed,values;color,linewidth,linestyle=up ? :solid : :dash,label=false)
                    end
                end
                for (values,marker) in ((d.lower,:dtriangle),(d.upper,:utriangle))
                    scatter!(p,d.speed,values;marker,markersize,markercolor=fill,
                        markerstrokecolor=color,markerstrokewidth=0.65,label=false)
                    bad=.!d.verified
                    any(bad) && scatter!(p,d.speed[bad],values[bad];marker=:xcross,
                        markersize=markersize+0.7,color=:black,markerstrokewidth=1.1,label=false)
                end
                if i == 1
                    label=(up ? "Upsweep" : "Downsweep")*" ($(round(case.dataset.angle;digits=2))°)"
                    plot!(p,[NaN],[NaN];color,linewidth,linestyle=up ? :solid : :dash,
                        marker=:circle,markersize,markercolor=fill,markerstrokecolor=color,label)
                end
            end
            if i == 1
                for (label,marker) in (("Minima",:dtriangle),("Maxima",:utriangle))
                    scatter!(p,[NaN],[NaN];marker,markersize,markercolor=:white,
                        markerstrokecolor=:black,markerstrokewidth=0.65,label)
                end
                any(c -> any(.!c.dataset.panels[i].verified),cases) && scatter!(p,[NaN],[NaN];
                    marker=:xcross,markersize=markersize+0.7,color=:black,markerstrokewidth=1.1,
                    label="Unverified sample")
            end
            isnothing(xlims) || Plots.xlims!(p,xlims)
            isnothing(ylims[i]) || Plots.ylims!(p,ylims[i])
            push!(panels,p)
            cap=plot(;xlims=(0,1),ylims=(0,1),axis=false,ticks=false,
                grid=false,framestyle=:none,margin=0Plots.mm,legend=false)
            annotate!(cap,0.5,0.5,text(Paper.PANEL_CAPTIONS[i],12,"Times Bold"))
            push!(captions,cap)
            if save_panels
                plot!(p;size=panel_size,dpi)
                for format in formats
                    path=joinpath(directory,"pazy_stepped_$(response).$format")
                    savefig(p,path)
                    push!(saved,path)
                end
            end
        end
        note=plot(;xlims=(0,1),ylims=(0,1),axis=false,ticks=false,grid=false,
            framestyle=:none,margin=0Plots.mm,legend=false)
        annotate!(note,0.5,0.5,text("Curves interpolate each sweep's final-window envelopes; they do not establish branch stability.",10,"Times"))
        figure=plot(panels[1],panels[2],captions[1],captions[2],panels[3],panels[4],captions[3],captions[4],note;
            layout=@layout([grid(4,2,heights=[0.44,0.06,0.44,0.06]); note{0.05h}]),
            size=figure_size,dpi,plot_title="Pazy stepped-speed bifurcation",plot_titlefontsize=14,fontfamily="Times")
        for format in formats
            path=joinpath(directory,"pazy_stepped_bifurcation.$format")
            savefig(figure,path)
            push!(saved,path)
        end
        csv=joinpath(directory,"interpolated_curves.csv")
        open(csv,"w") do io
            println(io,"direction,angle_of_attack_deg,response,airspeed_mps,minimum,maximum")
            for case in cases,d in case.dataset.panels,c in d.curves,k in eachindex(c.speed)
                println(io,join((case.direction,case.dataset.angle,d.response,c.speed[k],c.minimum[k],c.maximum[k]),','))
            end
        end
        push!(saved,csv)
        return saved
    finally
        Plots.backend(previous_backend)
    end
end
end # module

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("Usage: julia --project=. PazySteppedBifurcationPlots.jl path/to/stepped_run")
    ENV["GKSwstype"]=get(ENV,"GKSwstype","100")
    foreach(println,PazySteppedBifurcationPlots.plot_stepped_bifurcation(only(ARGS)))
end
