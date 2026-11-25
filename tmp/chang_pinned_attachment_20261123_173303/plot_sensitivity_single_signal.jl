include(joinpath(@__DIR__,"plot_wing_chordwise.jl"))

function save_sensitivity_curve(output,name,x,percent,xlabel,xlimits,xticks)
    query = sort!(unique!(vcat(collect(range(first(x),last(x);length=600)),x)))
    curve = interpolate_sensitivity(x,percent,query)
    all(isapprox.(interpolate_sensitivity(x,percent,x),percent;atol=1e-12)) ||
        error("Interpolation missed a sample in $name")
    minimum(percent)-1e-12 <= minimum(curve) <= maximum(curve) <= maximum(percent)+1e-12 ||
        error("Interpolation overshot a sample in $name")

    peak = 1.08max(maximum(percent),1e-3)
    tick = sensitivity_tick_step(peak)
    upper = tick * ceil(peak/tick)
    p = plot(query,curve;color=:black,linewidth=3,legend=false,grid=false,
        fontfamily="sans-serif",tickfontsize=10,labelfontsize=12,
        margin=5Plots.mm,framestyle=:axes,tick_direction=:in,
        xlabel,ylabel="Difference (%)",xlims=xlimits,xticks,ylims=(0,upper),
        yticks=0:tick:upper,size=(600,400),dpi=300)
    savefig(p,joinpath(output,name * ".pdf"))
    println(name," sampled differences: ",collect(zip(x,percent)))
end

function saved_mesh(directory,family,metrics,speed)
    data = Saved.load_current_results(directory;metrics)
    rows = filter(row -> row.family == family,data.rows)
    data.settings["physical"]["speed_mps"] == speed || error("Unexpected airspeed in $directory")
    [row.x for row in rows] == data.settings["sweeps"][family] || error("Incomplete $family sweep")
    all(row -> row.valid && all(isfinite,row.values),rows) || error("Invalid $family fit")
    return rows,data.settings["damping"]["minimum_fit_r_squared"]
end

function save_single_signal_figures(output,stepPath)
    results = joinpath(root,"lib","WingPropellerUVLM","examples","chang_linear_aeroelastic","output")
    spanDirectory = joinpath(results,"convergence_pinned_full_span_10threads_xaJCcF")
    chordDirectory = joinpath(results,"convergence_aeroelastic_10threads_20261123_163633_rDWuCC")
    propellerDirectory = joinpath(results,"convergence_aeroelastic_10threads_20261120_222350_aR9lDM")
    spanRows,_ = saved_mesh(spanDirectory,"wing_span",WING_METRICS,65.0)
    chordRows,minimumR2 = saved_mesh(chordDirectory,"wing_chord",WING_METRICS,65.0)
    chordRates = refit_common_window(chordDirectory,chordRows,"torsion",minimumR2)
    chordPercent = Saved.relative_error_percent.(abs.(diff(chordRates)),chordRates[2:end])
    save_sensitivity_curve(output,"wing_spanwise_torsion_65mps",[row.x for row in spanRows[2:end]],
        [row.successive_percent[2] for row in spanRows[2:end]],"Spanwise panels",(0,83),0:10:80)
    save_sensitivity_curve(output,"wing_chordwise_torsion_65mps",[row.x for row in chordRows[2:end]],
        chordPercent,"Chordwise panels",(0,83),0:10:80)

    for (family,name,xlabel) in (
        ("prop_radial","propeller_spanwise_pitch_80mps","Spanwise panels"),
        ("prop_chord","propeller_chordwise_pitch_80mps","Chordwise panels")
    )
        rows,_ = saved_mesh(propellerDirectory,family,Saved.METRICS,80.0)
        save_sensitivity_curve(output,name,[row.x for row in rows[2:end]],
            [row.successive_percent[1] for row in rows[2:end]],xlabel,(0,83),0:10:80)
    end

    data,header = readdlm(stepPath,'\t',Float64;header=true)
    column = findfirst(==("CT_mean"),vec(header))
    isnothing(column) && error("Missing CT_mean in the azimuthal-step results")
    order = sortperm(data[:,1])
    steps,thrust = data[order,1],data[order,column]
    allunique(steps) && all(isfinite,thrust) || error("Invalid azimuthal-step results")
    percent = Saved.relative_error_percent.(abs.(diff(thrust)),thrust[1:end-1])
    save_sensitivity_curve(output,"azimuthal_step_thrust_coefficient",steps[2:end],percent,
        "Azimuthal step (deg.)",(0,30),0:5:30)
end

function save_chordwise_doubling(output)
    directory = joinpath(root,"lib","WingPropellerUVLM","examples","chang_linear_aeroelastic","output",
        "convergence_aeroelastic_10threads_20261123_163633_rDWuCC")
    rows,minimumR2 = saved_mesh(directory,"wing_chord",WING_METRICS,65.0)
    rates = refit_common_window(directory,rows,"torsion",minimumR2)
    counts = [5,10,20,40,80]
    selected = findall(row -> row.x in counts,rows)
    [rows[i].x for i in selected] == counts || error("Missing doubling levels")
    selectedRates = rates[selected]
    percent = Saved.relative_error_percent.(abs.(diff(selectedRates)),selectedRates[2:end])
    save_sensitivity_curve(output,"wing_chordwise_torsion_doubling_from5_65mps",counts[2:end],percent,
        "Chordwise panels",(0,83),vcat(0,counts))
end

if first(ARGS) == "--wing-chord-doubling"
    save_chordwise_doubling(only(ARGS[2:end]))
else
    save_single_signal_figures(ARGS...)
end
