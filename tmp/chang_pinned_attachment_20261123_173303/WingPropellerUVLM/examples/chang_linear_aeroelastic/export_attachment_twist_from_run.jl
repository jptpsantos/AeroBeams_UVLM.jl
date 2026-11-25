# Recover the attachment twist from a still-live Chang simulation, after checking
# that its response matches the existing organized text history.
using DelimitedFiles

function export_matching_attachment_history(chang_run;
    source=joinpath(@__DIR__,"output","multi_prop","single_propeller_history.txt"),
    destination=joinpath(dirname(source),"single_propeller_history_with_attachment.txt"))
    lines=readlines(source)
    header=findfirst(line -> startswith(strip(line),"time_s\t"),lines)
    isnothing(header) && error("Missing tab-separated time_s header")
    names=split(strip(lines[header]),'\t')
    required=["time_s","propeller_1_pitch_deg","propeller_1_yaw_deg","tip_displacement_m","tip_twist_deg"]
    indices=[findfirst(==(name),names) for name in required]
    any(isnothing,indices) && error("Source lacks required columns")
    data=readdlm(IOBuffer(join(lines[header+1:end],"\n")),'\t',Float64)[:,indices]
    results=chang_run.results
    reference=(results.time,rad2deg.(results.propeller_pitch[1]),rad2deg.(results.propeller_yaw[1]),
        results.tip_displacement,rad2deg.(results.tip_twist))
    n=size(data,1)
    all(length(column)==n for column in reference) || error("Run and saved history have different lengths")
    for i in 1:5
        all(isapprox.(data[:,i],reference[i];rtol=1e-9,atol=1e-10)) ||
            error("Run does not match saved history in $(required[i]); refusing to combine different runs")
    end
    params=chang_run.model.parameters
    node=params.prop_attach_nodes[1]
    ndof=params.ndof
    node >= 1 && ndof >= 4 || error("Invalid attachment node or structural DOF count")
    history=chang_run.solution.displacement_history
    length(history)>=n || error("Run does not retain all displacement states")
    index=(node-2)*ndof+4
    attachment=node==1 ? zeros(n) : [rad2deg(history[i][index]) for i in 1:n]
    all(isfinite,attachment) || error("Attachment rotation contains non-finite values")
    ispath(destination) && error("Destination already exists: $destination")
    open(destination,"w") do io
        println(io,"Chang single-propeller response history with recovered attachment twist")
        println(io,"source = $(abspath(source))")
        println(io,"propeller_1_attachment_node = $node")
        println(io,"\ntime_s\tpropeller_1_pitch_deg\tpropeller_1_yaw_deg\tpropeller_1_attachment_wing_twist_deg\ttip_displacement_m\ttip_twist_deg")
        for i in 1:n
            println(io,join((data[i,1],data[i,2],data[i,3],attachment[i],data[i,4],data[i,5]),'\t'))
        end
    end
    println("Matching attachment history written to $destination")
    return destination
end
