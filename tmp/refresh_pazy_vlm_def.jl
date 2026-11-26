example = joinpath(@__DIR__,"..","lib","WingPropellerUVLM","examples","pazy_aerobeams_coupled")
include(joinpath(example,"run_pazy_s0_steady_pitch_range.jl"))
include(joinpath(example,"run_pazy_swept_steady_pitch_range.jl"))
LinearAlgebra.BLAS.set_num_threads(1)

function save_comparison(directory,header,data)
    mkpath(directory)
    open(joinpath(directory,"steady_pitch_range.csv"),"w") do io
        writedlm(io,header,',')
        writedlm(io,data,',')
    end
end

for (oldDirectory,newDirectory,options,name) in (
    ("s0_steady_pitch_range",S0_BEAM_ROD_DIRECTORY,S0_BEAM_ROD_OPTIONS,"S0 beam-rod model"),
    ("s0_inertia_rod_steady_pitch_range",S0_STUDY_DIRECTORY,S0_MODEL_OPTIONS,"S0 beam with tip-rod inertia"))
    data,header = readdlm(joinpath(example,"output",oldDirectory,"steady_pitch_range.csv"),',',Float64; header=true)
    @assert size(data,1) == (options.beamWingtipRod ? 183 : 39) && all(data[:,9] .== 1)
    for pitch in S0_PITCH_RANGE
        x0 = Float64[]
        for i in findall(==(pitch),data[:,1])
            native = solve_s0_native_point(data[i,2],deg2rad(pitch),x0,options)
            data[i,5],data[i,6],x0 = native.oop,native.twist,native.x
        end
        println("$name: checked VLM-def loads at $pitch degrees")
        flush(stdout)
    end
    save_comparison(newDirectory,header,data)
    plot_s0_steady_pitch_range(newDirectory; mesh=options.mesh,modelLabel=name)
end

data,header = readdlm(joinpath(example,"output","swept_steady_benchmark","steady_pitch_range.csv"),',',Float64; header=true)
@assert size(data,1) == 96 && all(data[:,10] .== 1)
for sweep in SWEPT_ANGLES_DEG, pitch in SWEPT_PITCH_RANGE
    x0 = Float64[]
    rows = findall(i -> isapprox(data[i,1],sweep) && data[i,2] == pitch,axes(data,1))
    @assert length(rows) == length(SWEPT_AIRSPEEDS)
    for i in rows
        native = solve_swept_native_point(data[i,3],deg2rad(data[i,18]),deg2rad(sweep),x0)
        data[i,6],data[i,7],data[i,19],data[i,21],x0 = native.oop,native.twist,native.strain,native.aoa,native.x
    end
    println("S$sweep: checked VLM-def loads at $pitch degrees")
    flush(stdout)
end
save_comparison(SWEPT_STUDY_DIRECTORY,header,data)
plot_pazy_swept_steady_pitch_range()
