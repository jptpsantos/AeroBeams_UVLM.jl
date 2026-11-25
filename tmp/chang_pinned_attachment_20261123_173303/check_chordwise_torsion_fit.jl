using DelimitedFiles,LinearAlgebra,Printf,Statistics

BLAS.set_num_threads(1)

function fit_damped_sinusoid(time,twist)
    tau = time .- first(time)
    baseline = [ones(length(tau)) tau]
    function residual(rate,frequency)
        envelope = exp.(rate .* tau)
        phase = 2frequency .* tau
        design = [baseline envelope .* cospi.(phase) envelope .* sinpi.(phase)]
        return sum(abs2,twist .- design * (design \ twist))
    end
    best = findmin([residual(rate,frequency)
        for rate in -2.5:0.125:0,frequency in 3.5:0.05:5.5])
    rate = (-2.5:0.125:0)[best[2][1]]
    frequency = (3.5:0.05:5.5)[best[2][2]]
    errorSum = best[1]
    for level in 1:6
        rates = range(rate-0.125/4^(level-1),rate+0.125/4^(level-1);length=9)
        frequencies = range(frequency-0.05/4^(level-1),frequency+0.05/4^(level-1);length=9)
        errorSum,index = findmin([residual(a,b) for a in rates,b in frequencies])
        rate,frequency = rates[index[1]],frequencies[index[2]]
    end
    rSquared = 1 - errorSum/sum(abs2,twist .- mean(twist))
    return rate,frequency,rSquared
end

directory = only(ARGS)
matrix,_ = readdlm(joinpath(directory,"case_matrix.csv"),',';header=true)
levels = Int.(matrix[:,5])
histories = map(string.(matrix[:,3])) do label
    data,header = readdlm(joinpath(directory,label,"response_history.csv"),',';header=true)
    columns = [findfirst(==(name),vec(header)) for name in ("time_s","tip_twist_deg")]
    return Float64.(data[:,Int.(columns)])
end
for (startTime,endTime) in ((0.87,2.78),(1.02,2.78),(0.87,2.60))
    println("Common free-response window ",startTime," through ",endTime," s")
    previous = NaN
    for (panels,data) in zip(levels,histories)
        mask = startTime .<= data[:,1] .<= endTime
        rate,frequency,rSquared = fit_damped_sinusoid(data[mask,1],data[mask,2])
        @printf("%2d panels  λ=% .8f  f=%.6f  R²=%.6f   difference=%.8f s⁻¹, %.6f%%\n",
            panels,rate,frequency,rSquared,abs(rate-previous),100abs(rate-previous)/abs(rate))
        previous = rate
    end
end
