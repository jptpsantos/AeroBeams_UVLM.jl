using Plots
using FFTW
using Polynomials
using Peaks
using Statistics

"""
    moving_block_damping(t, x; kwargs...)

Applies the moving block approach to estimate damping and dominant frequency
from a given time history.

Returns:
- `lambda_approx`: exponential growth/decay rate, in 1/s
- `frequency_approx`: dominant frequency, in Hz
- `omega_approx`: dominant angular frequency, in rad/s
"""
function moving_block_damping(t::AbstractVector, x::AbstractVector;
                              block_size::Int = 512,
                              size_ratio_lb::Float64 = 0.25,
                              size_ratio_ub::Float64 = 0.50,
                              peak_from_start::Int = 1,
                              peak_from_end::Int = 0)

    if length(t) != length(x)
        error("Time vector and signal vector must have the same length.")
    end

    peak_min = maximum(abs.(x)) / 100.0

    locs_all = argmaxima(x)
    locs = filter(i -> x[i] >= peak_min, locs_all)

    if length(locs) < (peak_from_start + peak_from_end + 1)
        error("Not enough peaks found in the signal. Try adjusting peak_from_start/end.")
    end

    start_idx = locs[peak_from_start]
    end_idx = locs[end - peak_from_end]

    t_test = t[start_idx:end_idx]
    x_test = x[start_idx:end_idx]

    sample_size = length(t_test)
    size_ratio = block_size / sample_size

    while size_ratio < size_ratio_lb
        block_size *= 2
        size_ratio = block_size / sample_size
    end

    while size_ratio > size_ratio_ub
        block_size = div(block_size, 2)
        size_ratio = block_size / sample_size
    end

    n_blocks = sample_size - block_size + 1

    if n_blocks < 1
        error("Block size is larger than the analyzed signal segment.")
    end

    dt = mean(diff(t_test))
    fs = 1.0 / dt

    lnXbar = zeros(Float64, n_blocks)
    freq_block = zeros(Float64, n_blocks)

    for i in 1:n_blocks
        x_i = x_test[i : i + block_size - 1]

        # Remove DC offset so the 0 Hz bin does not dominate.
        x_i = x_i .- mean(x_i)

        X_i = fft(x_i)
        Xbar_i = abs.(X_i ./ block_size)

        half_idx = div(block_size, 2) + 1
        Xbar_i = Xbar_i[1:half_idx]
        Xbar_i[2:end-1] .*= 2.0

        freqs = (0:half_idx-1) .* fs ./ block_size

        # Ignore DC bin when finding dominant oscillation.
        peak_idx = argmax(Xbar_i[2:end]) + 1

        lnXbar[i] = log(max(Xbar_i[peak_idx], eps()))
        freq_block[i] = freqs[peak_idx]
    end

    tau = t_test[1:n_blocks]

    p_fit = fit(tau, lnXbar, 1)
    coeff = coeffs(p_fit)

    lambda_approx = coeff[2]
    frequency_approx = median(freq_block)
    omega_approx = 2.0 * pi * frequency_approx

    p_test = plot(t_test, x_test,
                  linestyle=:dash,
                  linewidth=1,
                  label="Analyzed Segment",
                  title="Cut Signal",
                  xlabel="Time (s)")

    p_fit_plot = plot(tau, lnXbar,
                      linewidth=1,
                      label="lnXbar")

    plot!(p_fit_plot, tau, p_fit.(tau),
          linestyle=:dash,
          linewidth=1,
          label="Linear Fit",
          xlabel="Block start time (s)",
          ylabel="Log of moving block function (-)",
          title="Moving Block Fit")

    display(plot(p_test, p_fit_plot, layout=(2, 1), size=(800, 600)))

    return lambda_approx, frequency_approx, omega_approx
end

estimated_lambda, estimated_frequency, estimated_omega = moving_block_damping(t, t_prop)

println("Calculated lambda: ", estimated_lambda, " 1/s")
println("Calculated frequency: ", estimated_frequency, " Hz")
println("Calculated omega: ", estimated_omega, " rad/s")


"""
    find_flutter_speed(v1, zeta1, v2, zeta2)

Calculates the flutter speed where damping crosses zero using linear interpolation.
"""
function find_flutter_speed(v1::Real, zeta1::Real, v2::Real, zeta2::Real)
    if sign(zeta1) == sign(zeta2)
        @warn "Damping does not change sign between these points. The result is an extrapolation."
    end

    V_f = v1 - zeta1 * ((v2 - v1) / (zeta2 - zeta1))

    return V_f
end

damp = [
79	-0.011534265
80	0.030646505
]

V_1 = damp[1, 1]
zeta_1 = damp[1, 2]

V_2 = damp[2, 1]
zeta_2 = damp[2, 2]

V_flutter = find_flutter_speed(V_1, zeta_1, V_2, zeta_2)

println("Estimated Flutter Speed: ", round(V_flutter, digits=2), " m/s")