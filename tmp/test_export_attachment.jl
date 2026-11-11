using Test, DelimitedFiles
include(joinpath(pwd(),"lib","WingPropellerUVLM","examples","chang_linear_aeroelastic","export_attachment_twist_from_run.jl"))
d=readdlm(joinpath(pwd(),"tmp","hub_pitch_test","synthetic_history.txt"),'\t',Float64;skipstart=2)
source=tempname(joinpath(pwd(),"tmp"))
destination=tempname(joinpath(pwd(),"tmp"))
open(source,"w") do io
    println(io,"time_s\tpropeller_1_pitch_deg\tpropeller_1_yaw_deg\ttip_displacement_m\ttip_twist_deg")
    for i in axes(d,1)
        println(io,join(d[i,[1,2,3,5,6]],'\t'))
    end
end
run=(results=(time=d[:,1],propeller_pitch=[deg2rad.(d[:,2])],propeller_yaw=[deg2rad.(d[:,3])],tip_displacement=d[:,5],tip_twist=deg2rad.(d[:,6])),
    model=(parameters=(prop_attach_nodes=[2],ndof=6),),
    solution=(displacement_history=[[0.0,0,0,deg2rad(d[i,4]),0,0] for i in axes(d,1)],))
export_matching_attachment_history(run;source,destination)
recovered=readdlm(destination,'\t',Float64;skipstart=5)
@test recovered[:,4] ≈ d[:,4]
bad=merge(run,(results=merge(run.results,(tip_twist=run.results.tip_twist .+ 0.1,)),))
@test_throws ErrorException export_matching_attachment_history(bad;source,destination=tempname(joinpath(pwd(),"tmp")))
println("Attachment recovery and mismatch checks passed")
