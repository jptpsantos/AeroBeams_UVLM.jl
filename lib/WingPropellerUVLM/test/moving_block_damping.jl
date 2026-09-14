module MovingBlockDampingTests
using Test, FFTW, Statistics, LinearAlgebra
include(joinpath(@__DIR__,"..","examples","chang_linear_aeroelastic",
    "studies","convergence","MovingBlockDamping.jl"))
using .MovingBlockDamping

# Independent full-FFT transcription of moving_block_2.jl's numerical steps.
# QR regression replaces Polynomials.fit(degree=1); no plotting/demo execution.
function reference(t,x)
    threshold=maximum(abs.(x))/100
    peaks=filter(i -> x[i]>=threshold,[i for i in 2:length(x)-1
        if x[i]>x[i-1] && x[i]>=x[i+1]])
    cut=first(peaks):last(peaks)
    tx=t[cut]; xx=x[cut]; n=length(xx); b=512
    while b/n<.25; b*=2; end
    while b/n>.5; b=div(b,2); end
    count=n-b+1; dt=mean(diff(tx))
    logs=zeros(count); frequencies=zeros(count)
    for i in 1:count
        window=xx[i:i+b-1]; window=window.-mean(window)
        amplitudes=abs.(fft(window)./b)[1:div(b,2)+1]
        amplitudes[2:end-1].*=2
        index=argmax(amplitudes[2:end])+1
        logs[i]=log(max(amplitudes[index],eps()))
        frequencies[i]=(index-1)/(b*dt)
    end
    coefficients=hcat(ones(count),tx[1:count])\logs
    return coefficients[2],median(frequencies),b,count
end

@testset "Reference moving-block FFT method" begin
    t=collect(0.:.002:6.)
    for lambda in (-.15,.10), frequency in (4.1,4.4)
        x=exp.(lambda.*t).*sin.(2pi*frequency.*t)
        expected=reference(t,x)
        result=moving_block_metrics(t,x)
        m=result.summary
        @test m.available
        @test m.moving_block_lambda_per_s≈expected[1] atol=1e-10
        @test m.frequency_hz≈expected[2] atol=1e-10
        @test m.block_size==expected[3]
        @test m.block_count==expected[4]
        @test m.moving_block_lambda_per_s≈lambda atol=.01
        @test m.omega_radps≈2pi*m.frequency_hz
        @test m.frequency_resolution_hz≈1/(m.block_size*.002)
        @test all(isapprox.(diff(result.diagnostics.block_start_s),.002;atol=1e-12))
        @test all(isapprox.(moving_block_damping(t,x),(expected[1],expected[2],2pi*expected[2]);atol=1e-10))
    end
    outside=moving_block_metrics(t,exp.(-.15t).*sin.(2pi*8t);
        frequency_band_hz=(3.,6.),minimum_fit_r_squared=.8).summary
    @test outside.available
    @test !outside.valid_for_convergence
    @test outside.frequency_hz>6
    @test occursin("band",outside.reason)
    @test !moving_block_metrics(t,zeros(length(t))).summary.available
    @test_throws DimensionMismatch moving_block_damping(t,[1.,2.])
    @test_throws ArgumentError moving_block_damping(t,zeros(length(t)))
    @test_throws ArgumentError moving_block_metrics(t,copy(t);block_size=0)
    @test_throws ArgumentError moving_block_metrics(t,copy(t);size_ratio_lb=.7,size_ratio_ub=.5)
    @test_throws ArgumentError moving_block_metrics(t,copy(t);peak_from_start=0)
    @test_throws ArgumentError moving_block_metrics([0.,.1,.3],[1.,2.,3.])
    @test_throws ArgumentError moving_block_metrics([0.,.1,.1],[1.,2.,3.])
    @test_throws ArgumentError moving_block_metrics([0.,.1,.2],[1.,NaN,3.])
    duration=moving_block_metrics(t,exp.(-.15t).*sin.(2pi*4t);
        block_duration_s=1.,block_overlap=.9).summary
    @test duration.block_size==500
    @test duration.moving_block_lambda_per_s≈-.15 atol=.01
end
end
