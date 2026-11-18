ENV["GKSwstype"] = "100"
import WingPropellerUVLM
using Plots, Serialization, Printf

# Reuse the saved computed vertices when adjusting only presentation.
output = joinpath(@__DIR__,"s20_wake_equal_scale")
(;grid,frames,alpha,U,dt,nc,ns,R) = deserialize(joinpath(output,"computed_frames.jls"))
function lattice_lines(positions)
    x,y,z = Float64[],Float64[],Float64[]
    function add_line(points)
        for point in points
            push!(x,point[1]);push!(y,point[2]);push!(z,point[3])
        end
        push!(x,NaN);push!(y,NaN);push!(z,NaN)
    end
    for row in axes(positions,1)
        add_line(positions[row,:])
    end
    for col in axes(positions,2)
        add_line(positions[:,col])
    end
    return x,y,z
end
source = read(joinpath(@__DIR__,"render_s20_wake_equal_scale.jl"),String)
start = findfirst("    te_root = grid[:,end,1]",source)
stop = findfirst("    loads = UVLM.dimensional_loads(problem)",source)
plotcode = source[first(start):prevind(source,first(stop))]
plotcode = replace(plotcode,"left_margin=5Plots.mm"=>"left_margin=13Plots.mm")
include_string(@__MODULE__,plotcode,"s20_cached_plot.jl")
println("Replotted cached UVLM vertices with wider label margins.")
