# Analyze completed CSVs; this file never runs the aeroelastic solver.
using DelimitedFiles, LinearAlgebra, Statistics, Printf, FFTW
import Plots

function pazy_diagnostic_csv(path)
    data, header = readdlm(path, ',', Float64; header=true)
    return Dict(String(name) => data[:,i] for (i,name) in enumerate(vec(header)))
end

function pazy_diagnostic_fit(t, values)
    centered = t .- mean(t)
    slope = dot(centered, values .- mean(values)) / dot(centered, centered)
    fitted = mean(values) .+ slope .* centered
    total = sum(abs2, values .- mean(values))
    return slope, total > 0 ? 1-sum(abs2, values .- fitted)/total : NaN
end

function analyze_pazy_LCO_diagnostic(;
    output_directory=joinpath(@__DIR__, "output", "LCO_constant_speed_diagnostic"))
    names = ("baseline", "unforced", "half_dt", "double_wake")
    labels = ("Pulse baseline", "Unforced", "Half timestep", "Double wake")
    cases = []
    for (name,label) in zip(names,labels)
        path = joinpath(output_directory, "U45p0_$name")
        history = pazy_diagnostic_csv(joinpath(path, "history.csv"))
        windows = pazy_diagnostic_csv(joinpath(path, "windows.csv"))
        @assert all(v -> all(isfinite,v), values(history)) "Nonfinite history in $name"
        t, strain = history["time_s"], history["root_longitudinal_strain_ue"]
        dt = t[2]-t[1]
        @assert all(isapprox.(diff(t),dt; rtol=1e-10))
        @assert all(==(45.0),history["airspeed_mps"])
        @assert maximum(history["coupling_geometry_residual"]) <= 1e-4
        @assert maximum(history["coupling_load_residual"]) <= 1e-4
        reconstructed = 1e6 .* (history["root_axial_strain"] .+
            0.001125 .* history["root_curvature_2_per_m"])
        @assert maximum(abs.(reconstructed .- strain)) < 1e-8
        centers = (windows["window_start_s"] .+ windows["window_end_s"]) / 2
        rms = windows["detrended_rms_ue"]
        early = findall(x -> 0.75 <= x <= 2.5, centers)
        early_slope, early_r2 = pazy_diagnostic_fit(centers[early],log.(rms[early]))
        tail = length(centers)-3:length(centers)
        tail_slope, _ = pazy_diagnostic_fit(centers[tail],log.(rms[tail]))
        # A one-second Hann spectrum gives about 1 Hz bin spacing. This is a
        # finite-window dominant frequency, not a converged modal eigenvalue.
        last_second = findall(x -> x >= last(t)-1.0, t)
        local_time = t[last_second] .- mean(t[last_second])
        signal = strain[last_second] .- mean(strain[last_second])
        signal .-= local_time .* (dot(local_time,signal)/dot(local_time,local_time))
        hann = 0.5 .- 0.5 .* cos.(2pi .* (0:length(signal)-1) ./ (length(signal)-1))
        power = abs2.(rfft(signal .* hann))
        frequencies = (0:length(power)-1) ./ (length(signal)*dt)
        band = findall(f -> 1.0 <= f <= 100.0,frequencies)
        peak = band[argmax(power[band])]
        tail_rms = rms[end-2:end]
        tail_spread = 100 * (maximum(tail_rms)-minimum(tail_rms))/mean(tail_rms)
        push!(cases,(; name,label,history,windows,centers,dt,early_slope,early_r2,
            tail_slope,frequency=frequencies[peak],tail_spread))
    end
    baseline, unforced, half_dt, double_wake = cases
    @assert baseline.history["root_longitudinal_strain_ue"][1] ==
        unforced.history["root_longitudinal_strain_ue"][1]
    @assert all(iszero,unforced.history["tip_force_A1_N"])
    @assert maximum(abs.(baseline.history["tip_force_A1_N"] .-
        double_wake.history["tip_force_A1_N"])) < 1e-12
    @assert maximum(abs.(baseline.history["tip_force_A1_N"] .-
        half_dt.history["tip_force_A1_N"][1:2:end])) < 1e-12

    metrics_header = ["case", "dt_s", "initial_strain_ue", "initial_tip_OOP_m",
        "final_window_half_range_ue", "final_window_RMS_ue", "early_log_RMS_growth_per_s",
        "early_log_RMS_fit_R2", "tail_log_RMS_growth_per_s",
        "last_three_window_RMS_spread_pct", "last_second_dominant_frequency_Hz"]
    rows = [Any[c.name,c.dt,first(c.history["root_longitudinal_strain_ue"]),
        first(c.history["tip_out_of_plane_m"]),last(c.windows["detrended_half_range_ue"]),
        last(c.windows["detrended_rms_ue"]),c.early_slope,c.early_r2,c.tail_slope,
        c.tail_spread,c.frequency] for c in cases]
    writedlm(joinpath(output_directory,"analysis_metrics.csv"),
        vcat(permutedims(metrics_header),permutedims(hcat(rows...))),',')

    colors = (:dodgerblue,:orangered,:forestgreen,:purple)
    growth = Plots.plot(; xlabel="Time [s]",ylabel="Window RMS strain [microstrain]",
        yscale=:log10,legend=:topleft,title="Post-pulse growth (0.25 s windows)",
        margin=7Plots.mm)
    for (c,color) in zip(cases,colors)
        Plots.plot!(growth,c.centers,c.windows["detrended_rms_ue"];
            color,label=c.label,linewidth=1.8,marker=:circle,markersize=3)
    end
    equilibria = Plots.bar(1:4,[first(c.history["root_longitudinal_strain_ue"]) for c in cases];
        xticks=(1:4,["Pulse\nbaseline","Unforced","Half\ntimestep","Double\nwake"]),
        color=collect(colors),legend=false,ylabel="Initial surface strain [microstrain]",
        title="Equilibrium changes with timestep",margin=7Plots.mm)
    comparison = Plots.plot(growth,equilibria; layout=(1,2),size=(1350,550),dpi=180)
    Plots.savefig(comparison,joinpath(output_directory,"analysis_comparison.png"))
    traces = map(cases,colors) do c,color
        indices = findall(t -> t >= last(c.history["time_s"])-0.25,c.history["time_s"])
        values = c.history["root_longitudinal_strain_ue"][indices]
        Plots.plot(c.history["time_s"][indices],values .- mean(values);
            color,label=false,xlabel="Time [s]",ylabel="Strain [microstrain]",
            title="$(c.label): final 0.25 s, mean removed",titlefontsize=10,
            margin=5Plots.mm)
    end
    Plots.savefig(Plots.plot(traces...; layout=(2,2),size=(1200,750),dpi=180),
        joinpath(output_directory,"analysis_late_responses.png"))

    report_path = joinpath(output_directory,"analysis_report.md")
    open(report_path,"w") do io
        println(io,"# Pazy constant-speed diagnostic: completed results\n")
        println(io,"All four cases completed at 45 m/s, 5 degrees root incidence, about 5 s physical time, with the same 10 x 20 surface mesh. No new simulations were run for this analysis.\n")
        println(io,"**Main result: the present response is strongly timestep-sensitive. The small downsweep cannot be attributed solely to missing excitation, and the current results do not establish numerical convergence.**\n")
        println(io,"| Case | dt [ms] | Initial strain [microstrain] | Final-window half range [microstrain] | Early RMS growth [1/s] | Tail RMS growth [1/s] | Late dominant frequency [Hz] |")
        println(io,"|---|---:|---:|---:|---:|---:|---:|")
        for c in cases
            @printf(io,"| %s | %.6f | %.2f | %.3f | %.3f | %.3f | %.2f |\n",
                c.label,1e3*c.dt,first(c.history["root_longitudinal_strain_ue"]),
                last(c.windows["detrended_half_range_ue"]),c.early_slope,c.tail_slope,c.frequency)
        end
        println(io,"\nHalf range and RMS come from a linear-detrended 0.25 s window; the final complete window is 4.64-4.89 s. The early growth fit uses window centers between 0.75 and 2.5 s; it is a finite-response estimate, not an eigenvalue. Tail growth uses the last four windows and may include the approach to saturation. Frequencies are Hann FFT peaks over the last 1 s, with about 1 Hz frequency-bin spacing.\n")
        println(io,"## Interpretation\n")
        println(io,"1. **Excitation and finite time matter.** The unforced response grows from very small disturbances to about 109 microstrain by the final window. The pulse case reaches about 1357 microstrain over the same elapsed time. Both grow with the baseline timestep; the unforced case has not reached a plateau. Residual equilibrium mismatch/roundoff supplies disturbances in the nominally unforced run. This supports delayed growth as a contributor to the ramp result, but does not quantify its contribution to the downsweep.\n")
        println(io,"2. **Timestep sensitivity dominates this comparison.** Halving dt reduces the final-window half range by about 99.92%, while reducing the magnitude of the initial strain by about 10.21%. Initial tip displacement changes from about 0.11990 m to 0.10818 m. The refined response is still growing slowly (about 0.29 1/s), so 1.06 microstrain is a transient value, not its settled LCO amplitude. Because the equilibrium changes too, the dynamic runs are linearizing/evolving about different configurations.\n")
        println(io,"3. **A longer wake has a smaller effect here, but wake convergence is unproven.** Doubling retained wake length from 6.49 to 12.98 chords changes the final-window half range by about -3.51% and initial strain magnitude by about +1.25%. It also changes the growth trajectory. The double-wake case is nearly stationary over its final four RMS windows; this short interval alone does not establish a unique, fully settled LCO.\n")
        println(io,"4. **Iteration convergence passed; discretization convergence did not.** Recorded geometry and load residuals remain within 1e-4, and the static initializations report convergence. Those tests measure iteration consistency at each chosen discretization. They do not remove the differences between timestep/wake discretizations.\n")
        println(io,"## What this means for the sweeps\n")
        println(io,"The old downsweep result remains a computed transient for its selected settings. These tests do not establish a plotting bug, a unique solver bug, or agreement with the reference's LCO branch. Raising the downsweep starting speed also reduces the shared timestep in the current driver, so a comparison against the old 52 m/s case changes both starting speed and numerical resolution. To isolate starting-speed effects, rerun both configurations at the same dt and physical wake length after identifying suitable discretization settings.\n")
        println(io,"## Recommended order of further work\n")
        println(io,"1. First investigate fixed-geometry aerodynamic/steady-load sensitivity to dt, without structural motion or a pulse. Check circulation and total loads together with wake geometry. The roughly 10% equilibrium-strain shift already exists before the pulse.\n")
        println(io,"2. Inspect near-trailing-edge wake spacing and finite-core treatment. In the current code the shedding offset is eta*V*dt. At eta=0.1 and 45 m/s this is about 2.14 mm for baseline and 1.07 mm at half dt, compared with the fixed 1 mm core radius. This is a concrete numerical sensitivity to investigate, not a confirmed cause.\n")
        println(io,"3. Then use an additional timestep refinement and controlled wake/core/mesh checks to establish trends in steady loads and small-disturbance growth. Keep physical wake length and pulse fixed; tighter coupling/Newton tolerances may also need checking for small-amplitude responses.\n")
        println(io,"4. Once those quantities are sufficiently insensitive, extend the constant-speed runs as needed for amplitude stationarity and repeat speed ramps with controlled excitation/initial conditions. A 45 m/s check alone cannot identify the full flutter boundary or reproduce the reference diagram.\n")
        println(io,"## Files and checks\n")
        println(io,"[Comparison figure](analysis_comparison.png), [late response traces](analysis_late_responses.png), [computed metrics](analysis_metrics.csv), [original summary](diagnostic_summary.csv).\n")
        println(io,"The analysis checked finite CSV values, uniform time grids, constant 45 m/s speed, reconstructed surface strain, and recorded coupling residuals. Baseline/unforced initial strains match exactly. Forced cases have matching pulse values at common sample times. Baseline and half dt retain the same nominal physical wake age/length (30 and 60 rows respectively); double wake uses 60 rows at baseline dt.\n")
    end
    println("Analysis written to $report_path")
    for c in cases
        @printf("%s: early growth %.5f /s (R2 %.5f), tail growth %.5f /s, late frequency %.3f Hz, last-three RMS spread %.3f%%\n",
            c.name,c.early_slope,c.early_r2,c.tail_slope,c.frequency,c.tail_spread)
    end
    return (; cases,report_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    analyze_pazy_LCO_diagnostic()
end
