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
function moving_block_damping(t::AbstractVector,x::AbstractVector;
    block_size::Int=512,
    size_ratio_lb::Float64=0.25,
    size_ratio_ub::Float64=0.50,
    peak_from_start::Int=1,
    peak_from_end::Int=0)

    if length(t) != length(x)
        error("Time vector and signal vector must have the same length.")
    end

    peakMin = maximum(abs.(x)) / 100.0

    locsAll = argmaxima(x)
    locs = filter(i -> x[i] >= peakMin,locsAll)

    if length(locs) < (peak_from_start + peak_from_end + 1)
        error("Not enough peaks found in the signal. Try adjusting peak_from_start/end.")
    end

    startIdx = locs[peak_from_start]
    endIdx = locs[end-peak_from_end]

    tTest = t[startIdx:endIdx]
    xTest = x[startIdx:endIdx]

    sampleSize = length(tTest)
    sizeRatio = block_size / sampleSize

    while sizeRatio < size_ratio_lb
        block_size *= 2
        sizeRatio = block_size / sampleSize
    end

    while sizeRatio > size_ratio_ub
        block_size = div(block_size,2)
        sizeRatio = block_size / sampleSize
    end

    nBlocks = sampleSize - block_size + 1

    if nBlocks < 1
        error("Block size is larger than the analyzed signal segment.")
    end

    dt = mean(diff(tTest))
    fs = 1.0 / dt

    lnXbar = zeros(Float64,nBlocks)
    freqBlock = zeros(Float64,nBlocks)

    for i = 1:nBlocks
        x_i = xTest[i:i+block_size-1]

        # Remove DC offset so the 0 Hz bin does not dominate.
        x_i = x_i .- mean(x_i)

        X_i = fft(x_i)
        Xbar_i = abs.(X_i ./ block_size)

        halfIdx = div(block_size,2) + 1
        Xbar_i = Xbar_i[1:halfIdx]
        Xbar_i[2:end-1] .*= 2.0

        freqs = (0:halfIdx-1) .* fs ./ block_size

        # Ignore DC bin when finding dominant oscillation.
        peakIdx = argmax(Xbar_i[2:end]) + 1

        lnXbar[i] = log(max(Xbar_i[peakIdx],eps()))
        freqBlock[i] = freqs[peakIdx]
    end

    tau = tTest[1:nBlocks]

    p_fit = fit(tau,lnXbar,1)
    coeff = coeffs(p_fit)

    lambda_approx = coeff[2]
    frequency_approx = median(freqBlock)
    omega_approx = 2.0 * pi * frequency_approx

    pTest = plot(tTest,xTest,
        linestyle=:dash,
        linewidth=1,
        label="Analyzed Segment",
        title="Cut Signal",
        xlabel="Time (s)")

    pFitPlot = plot(tau,lnXbar,
        linewidth=1,
        label="lnXbar")

    plot!(pFitPlot,tau,p_fit.(tau),
        linestyle=:dash,
        linewidth=1,
        label="Linear Fit",
        xlabel="Block start time (s)",
        ylabel="Log of moving block function (-)",
        title="Moving Block Fit")

    display(plot(pTest,pFitPlot,layout=(2,1),size=(800,600)))

    return lambda_approx,frequency_approx,omega_approx
end

estimated_lambda,estimated_frequency,estimated_omega = moving_block_damping(t,tProp)

println("Calculated lambda: ",estimated_lambda," 1/s")
println("Calculated frequency: ",estimated_frequency," Hz")
println("Calculated omega: ",estimated_omega," rad/s")

########################################################################################

using Statistics

results = chang_run.results
config = chang_run.config

t = results.time

# End of the applied pitch impulse.
impulse_end = config.excitation.impulse_start_s +
              config.excitation.impulse_duration_s

# Initial fitting-window choice; adjust after inspecting the response.
fit_start = max(0.75,impulse_end + 0.10)
fit_end = last(t)

keep = (t .>= fit_start) .& (t .<= fit_end)
any(keep) || error("No samples in the fitting interval.")

t_fit = t[keep]

pitch_fit = results.propeller_pitch[1][keep]
yaw_fit = results.propeller_yaw[1][keep]

# Center before the positive-peak selection in your function.
pitch_fit .-= mean(pitch_fit)
yaw_fit .-= mean(yaw_fit)

mb_options = (
    block_size=4096,
    size_ratio_lb=0.25,
    size_ratio_ub=0.50,
    peak_from_start=1,
    peak_from_end=0
)

lambda_pitch,frequency_pitch,omega_pitch =
    moving_block_damping(t_fit,pitch_fit; mb_options...)

lambda_yaw,frequency_yaw,omega_yaw =
    moving_block_damping(t_fit,yaw_fit; mb_options...)

# Equivalent damping ratios, assuming one dominant exponential oscillation.
zeta_pitch = -lambda_pitch / hypot(lambda_pitch,omega_pitch)
zeta_yaw = -lambda_yaw / hypot(lambda_yaw,omega_yaw)

println("Pitch: lambda = $lambda_pitch 1/s, f = $frequency_pitch Hz")
println("Yaw:   lambda = $lambda_yaw 1/s, f = $frequency_yaw Hz")
println("Pitch damping ratio = $(100*zeta_pitch)%")
println("Yaw damping ratio   = $(100*zeta_yaw)%")

############################

"""
    find_flutter_speed(v1, zeta1, v2, zeta2)

Calculates the flutter speed where damping crosses zero using linear interpolation.
"""
function find_flutter_speed(v1::Real,zeta1::Real,v2::Real,zeta2::Real)
    if sign(zeta1) == sign(zeta2)
        @warn "Damping does not change sign between these points. The result is an extrapolation."
    end

    V_f = v1 - zeta1 * ((v2 - v1) / (zeta2 - zeta1))

    return V_f
end

damp = [
    64 -0.52260291
    65 0.020612797
]

V_1 = damp[1,1]
zeta_1 = damp[1,2]

V_2 = damp[2,1]
zeta_2 = damp[2,2]

V_flutter = find_flutter_speed(V_1,zeta_1,V_2,zeta_2)

println("Estimated Flutter Speed: ",round(V_flutter,digits=2)," m/s")