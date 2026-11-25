# Response plotting helpers for the Chang example.

"""
    plot_chang_time_histories(...)

Create the standard Chang response plot. A bare-wing case shows tip displacement
and twist in two panels. A single-propeller case shows the
wing-tip displacement and twist plus propeller pitch and yaw. A two-propeller
case shows the inboard and outboard pitch/yaw histories, with the spanwise
roles inferred from attachment eta.
"""
function plot_chang_time_histories(
    time,
    tip_displacement,
    tip_twist,
    propeller_pitch,
    propeller_yaw,
    propeller_eta;
    span_length::Real,
    time_limit_s::Real,
    output_path::AbstractString
)
    number_of_propellers = length(propeller_pitch)
    length(propeller_yaw) == number_of_propellers ||
        throw(DimensionMismatch("Pitch and yaw histories must have the same propeller count"))
    length(propeller_eta) == number_of_propellers ||
        throw(DimensionMismatch("One attachment eta is required per propeller"))
    time_limit_s > 0 || throw(ArgumentError("Plot time limit must be positive"))

    common_style = (
        linewidth=6,
        color=:black,
        tickfontsize=20,
        labelfontsize=24,
        left_margin=22Plots.mm,
        bottom_margin=12Plots.mm,
        top_margin=5Plots.mm,
        right_margin=5Plots.mm,
        grid=false,
        label=false
    )
    timeLimits = (0.0,Float64(time_limit_s))

    if number_of_propellers == 0
        wing_style = merge(common_style,(;
            right_margin=15Plots.mm,
            top_margin=10Plots.mm
        ))
        panels = [
            Plots.plot(time,-tip_displacement ./ span_length .* 100;
                ylabel="Tip displacement (%)",xlims=timeLimits,wing_style...),
            Plots.plot(time,rad2deg.(tip_twist);
                xlabel="Time (s)",ylabel="Tip twist angle (deg)",
                xlims=timeLimits,wing_style...)
        ]
    elseif number_of_propellers == 1
        panels = [
            Plots.plot(time,-tip_displacement ./ span_length .* 100;
                ylabel="Tip displacement (%)",xlims=timeLimits,
                ylims=(-0.3,0.3),yticks=-0.3:0.1:0.3,common_style...),
            Plots.plot(time,rad2deg.(tip_twist);
                ylabel="Tip twist angle (deg)",xlims=timeLimits,
                ylims=(-0.2,0.2),common_style...),
            Plots.plot(time,rad2deg.(propeller_pitch[1]);
                ylabel="Propeller pitch (deg)",xlims=timeLimits,
                ylims=(-7.5,7.5),common_style...),
            Plots.plot(time,rad2deg.(propeller_yaw[1]);
                xlabel="Time (s)",ylabel="Propeller yaw (deg)",
                xlims=timeLimits,ylims=(-7.5,7.5),common_style...)
        ]
    elseif number_of_propellers == 2
        inboard_index,outboard_index = sortperm(propeller_eta)
        println(
            "Plotting inboard P$(inboard_index) (eta=$(propeller_eta[inboard_index])) " *
            "and outboard P$(outboard_index) (eta=$(propeller_eta[outboard_index])) responses..."
        )
        panels = [
            Plots.plot(time,rad2deg.(propeller_pitch[inboard_index]);
                ylabel="Inboard pitch (deg)",xlims=timeLimits,
                ylims=(-7.5,7.5),common_style...),
            Plots.plot(time,rad2deg.(propeller_yaw[inboard_index]);
                ylabel="Inboard yaw (deg)",xlims=timeLimits,
                ylims=(-7.5,7.5),common_style...),
            Plots.plot(time,rad2deg.(propeller_pitch[outboard_index]);
                ylabel="Outboard pitch (deg)",xlims=timeLimits,
                ylims=(-7.5,7.5),common_style...),
            Plots.plot(time,rad2deg.(propeller_yaw[outboard_index]);
                xlabel="Time (s)",ylabel="Outboard yaw (deg)",
                xlims=timeLimits,ylims=(-7.5,7.5),common_style...)
        ]
    else
        @warn "Automatic Chang time-history plotting supports zero, one, or two propellers" number_of_propellers
        return nothing
    end

    finalPlot = Plots.plot(panels...; layout=(length(panels),1),
        size=(1000,500 * length(panels)))
    Plots.savefig(finalPlot,output_path)
    display(finalPlot)
    println("Time-history plot written to $output_path")
    return finalPlot
end
