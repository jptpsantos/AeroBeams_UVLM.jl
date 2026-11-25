using Printf,Statistics

root = normpath(joinpath(@__DIR__,"..",".."))
convergence = joinpath(root,"lib","WingPropellerUVLM","examples",
    "chang_linear_aeroelastic","studies","convergence")
include(joinpath(convergence,"plot_chang_convergence_sensitivity.jl"))

function refit_common_window(directory,rows,signal,minimumR2)
    blocks = map(rows) do row
        file = joinpath(directory,row.label,"wing_$(signal)_moving_block.csv")
        Saved.read_table(file)
    end
    times = [[Saved.number(block["block_start_s"]) for block in series] for series in blocks]
    startTime = maximum(first,times)
    endTime = minimum(last,times)
    indices = [findall(t -> startTime <= t <= endTime,time) for time in times]
    commonTime = first(times)[first(indices)]
    length(commonTime) >= 2 || error("Insufficient shared moving blocks for $signal")
    centeredTime = commonTime .- mean(commonTime)
    @printf("%s common block starts: %.9f to %.9f s (%d blocks)\n",
        signal,startTime,endTime,length(commonTime))

    return map(eachindex(rows)) do i
        times[i][indices[i]] == commonTime || error("Different block times for $(rows[i].label)")
        logAmplitude = [Saved.number(block["log_amplitude"]) for block in blocks[i][indices[i]]]
        centeredAmplitude = logAmplitude .- mean(logAmplitude)
        lambda = sum(centeredTime .* centeredAmplitude) / sum(abs2,centeredTime)
        rSquared = 1 - sum(abs2,centeredAmplitude .- lambda .* centeredTime) / sum(abs2,centeredAmplitude)
        rSquared >= minimumR2 || error("Common-window fit below the R² threshold for $(rows[i].label)")
        @printf("%s N=%d: lambda %.9f, R2 %.6f\n",signal,rows[i].x,lambda,rSquared)
        return lambda
    end
end

function plot_wing_chordwise(args)
    directory,figurePath = args[1:2]
    data = Saved.load_current_results(directory;metrics=WING_METRICS)
    rows = filter(r -> r.family == "wing_chord",data.rows)
    @assert [r.x for r in rows] == data.settings["sweeps"]["wing_chord"]
    all(r -> r.valid && all(isfinite,r.values),rows) || error("Incomplete or invalid damping fits")

    if args[3:end] == ["--common-window"]
        minimumR2 = data.settings["damping"]["minimum_fit_r_squared"]
        bending,torsion = map(signal -> refit_common_window(directory,rows,signal,minimumR2),
            ("bending","torsion"))
        rates = collect(zip(bending,torsion))
        rows = map(eachindex(rows)) do i
            previous = i == 1 ? (NaN,NaN) : abs.(rates[i] .- rates[i-1])
            successivePercent = Saved.relative_error_percent.(previous,rates[i])
            return (; x=rows[i].x,previous,successive_percent=successivePercent)
        end
    end

    for (coarse,fine) in zip(rows[1:end-1],rows[2:end])
        @printf("%2d to %2d | damping difference %.8f, %.8f | percent %.5f, %.5f\n",
            coarse.x,fine.x,fine.previous[1:2]...,fine.successive_percent[1:2]...)
    end

    p = sensitivity_panel(rows,"wing_chord","Chordwise panels",:percent,WING_SIGNALS;showLegend=false)
    plot!(p;ylabel="Difference (%)",xlims=(0,last(rows).x+3))
    savefig(p,figurePath)
end

if abspath(PROGRAM_FILE) == @__FILE__
    plot_wing_chordwise(ARGS)
end
