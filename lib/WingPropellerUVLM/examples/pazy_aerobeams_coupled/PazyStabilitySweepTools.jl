# Small helpers for CSV output and the documented comparison plot.
using DelimitedFiles, Plots

sweep_number(x) = replace(replace(string(Float64(x)), r"\.0$" => ""), "." => "p", "-" => "m")

function write_sweep_csv(path, header, data)
    open(path,"w") do io
        println(io,join(header,','))
        writedlm(io,data,',')
    end
end

# Positive growth means instability. NaN breaks a bracket: never skip a failed case.
function sweep_crossings(speeds, growth)
    crossings = NamedTuple[]
    for j in 2:length(speeds)
        g1, g2 = growth[j-1], growth[j]
        isfinite(g1) && isfinite(g2) && g1*g2 < 0 || continue
        speed = speeds[j-1] - g1*(speeds[j]-speeds[j-1])/(g2-g1)
        push!(crossings,(;speed,kind=g1 < 0 ? "onset" : "offset"))
    end
    return crossings
end

# Cached TRACKED mode 3, as used by PazyWingFlutterPitchRange.jl.
# This reads the existing eigenanalysis, not a new strip-theory time-domain run.
function sweep_aerobeams_boundary()
    root = normpath(joinpath(@__DIR__,"..","..","..",".."))
    angles = unique(vcat(collect(0:0.25:1),collect(1:0.5:7)))
    onset, offset = fill(NaN,length(angles)), fill(NaN,length(angles))
    for (i,angle) in enumerate(angles)
        file = joinpath(root,"test","newTestDataGenerators","PazyWingFlutterPitchRange","damps$i.txt")
        damping = readdlm(file,Float64)
        for crossing in sweep_crossings(collect(0:120),damping[:,3])
            if crossing.kind == "onset" && isnan(onset[i])
                onset[i] = crossing.speed
            elseif crossing.kind == "offset" && isnan(offset[i])
                offset[i] = crossing.speed
            end
        end
    end
    return angles,onset,offset
end

function plot_sweep_comparison(boundaries, sweep_angles;
    aerobeams_color, literature_colors, labels, bracket_rows=Any[],
    aerobeams_styles=(:solid,:dot),
    experimental_colors=(:red,:green), uvlm_color=:red,
    uvlm_fill_alpha=0.25, line_width=2, marker_size=10,
    literature_marker_size=4, marker_stroke_width=0)
    gr()
    root = normpath(joinpath(@__DIR__,"..","..","..",".."))
    angles,aerobeams_onset,aerobeams_offset = sweep_aerobeams_boundary()
    sweep_angles = sort(unique(sweep_angles))
    p = plot()

    # For a bracketed offset, shade only up to the last confirmed growing speed.
    uvlm_onset = fill(NaN,length(sweep_angles))
    uvlm_offset = fill(NaN,length(sweep_angles))
    for (i,angle) in enumerate(sweep_angles)
        onsets = sort([row[2] for row in boundaries if row[1] == angle && row[3] == "onset"])
        offsets = sort([row[2] for row in boundaries if row[1] == angle && row[3] == "offset"])
        isempty(onsets) && continue
        offset_index = findfirst(speed -> speed > first(onsets),offsets)
        bracket_offsets = sort([row[2] for row in bracket_rows if
            row[1] == angle && row[4] == "offset" && row[2] > first(onsets)])
        offset = isnothing(offset_index) ?
            (isempty(bracket_offsets) ? NaN : first(bracket_offsets)) : offsets[offset_index]
        uvlm_onset[i],uvlm_offset[i] = first(onsets),offset
    end
    # Use the same single-polygon shape as the AeroBeams hump region. A gap
    # between valid pitch angles starts a new polygon instead of being bridged.
    region_labeled = false
    run_start = nothing
    for i in 1:(length(sweep_angles)+1)
        available = i <= length(sweep_angles) && isfinite(uvlm_onset[i]) && isfinite(uvlm_offset[i])
        if available && isnothing(run_start)
            run_start = i
        elseif !available && !isnothing(run_start)
            run_end = i - 1
            if run_end > run_start
                indices = run_start:run_end
                plot!(p,Shape(vcat(uvlm_onset[indices],reverse(uvlm_offset[indices])),
                    vcat(sweep_angles[indices],reverse(sweep_angles[indices])));
                    fillcolor=plot_color(uvlm_color,uvlm_fill_alpha),
                    linecolor=:black,linewidth=line_width,
                    label=region_labeled ? false : labels.uvlm_region)
                region_labeled = true
            end
            run_start = nothing
        end
    end

    # Experimental points and literature curves from the AeroBeams Pazy plot.
    for (speeds,pitch,label,color,marker) in (
        ([49,43,38],[3,5,7],labels.test_onset_up,experimental_colors[1],:rtriangle),
        ([58,51,46],[3,5,7],labels.test_offset_up,experimental_colors[2],:rtriangle),
        ([55,48],[3,5],labels.test_onset_down,experimental_colors[1],:ltriangle),
        ([40,36],[3,5],labels.test_offset_down,experimental_colors[2],:ltriangle))
        scatter!(p,speeds,pitch;label,color,marker,
            markersize=marker_size,markerstrokewidth=marker_stroke_width)
    end
    for (file,label,color,style,marker) in (
        ("Pazy/flutterBoundaryPitch_UMNAST.txt",labels.umnast_loss,literature_colors[1],:dash,:circle),
        ("Pazy/flutterBoundaryPitch_UMNAST_PanelCoeffs.txt",labels.umnast_panel,literature_colors[2],:dash,:circle),
        ("sweptPazy/flutterBoundary_UVsPitch_Lambda0_Sharpy.txt",labels.sharpy,literature_colors[3],:dash,:diamond))
        data = readdlm(joinpath(root,"test","referenceData",file),Float64)
        file == "sweptPazy/flutterBoundary_UVsPitch_Lambda0_Sharpy.txt" && (data = data[:,1:19])
        plot!(p,data[1,:],data[2,:];label,color,linestyle=style,
            linewidth=line_width,marker,
            markersize=literature_marker_size,markerstrokewidth=marker_stroke_width)
    end
    # Earlier AeroBeams eigenvalue results use one color and two line styles.
    plot!(p,aerobeams_onset,angles;label=labels.aerobeams_onset,color=aerobeams_color,
        linestyle=aerobeams_styles[1],linewidth=line_width)
    plot!(p,aerobeams_offset,angles;label=labels.aerobeams_offset,color=aerobeams_color,
        linestyle=aerobeams_styles[2],linewidth=line_width)
    return p
end
