# Chang wing–propeller aeroelastic response.
# Edit chang_case.jl for the case, numerical controls, and output.
# CHANG_* environment variables override those defaults when config is loaded.

import Pkg
let project = normpath(joinpath(@__DIR__,"..",".."))
    if Base.active_project() != joinpath(project,"Project.toml")
        Pkg.activate(project)
    end
end

if !isdefined(@__MODULE__,:ChangAeroelastic) ||
   !isdefined(ChangAeroelastic,:load_chang_configuration)
    include(joinpath(@__DIR__,"src","ChangAeroelastic.jl"))
end
using .ChangAeroelastic

function run_chang_example()
    config = load_chang_configuration()
    return run_chang(config)
end

"""Write organized histories; `df` is a NamedTuple of columns in CSV order."""
function export_chang_organized_history(chang_run)
    config = chang_run.config
    results = chang_run.results
    eta = config.propeller.attachment_eta
    attachment_nodes = results.propeller_attachment_node

    if length(eta) == 1
        df = (;
            time_s=results.time,
            tip_displacement_m=results.tip_displacement,
            tip_twist_deg=rad2deg.(results.tip_twist),
            propeller_1_pitch_deg=rad2deg.(results.propeller_pitch[1]),
            propeller_1_yaw_deg=rad2deg.(results.propeller_yaw[1]),
            propeller_1_attachment_node=fill(attachment_nodes[1],length(results.time)),
            propeller_1_attachment_wing_twist_deg=rad2deg.(results.propeller_attachment_twist[1])
        )
        csvName = "single_propeller_history.csv"
        txtName = "single_propeller_history.txt"
        header = "Chang single-propeller response history"
        printlnText = "time_s\tpropeller_1_pitch_deg\tpropeller_1_yaw_deg\tpropeller_1_attachment_node\tpropeller_1_attachment_wing_twist_deg\ttip_displacement_m\ttip_twist_deg"
    elseif length(eta) == 2
        perm = sortperm(eta)
        inboard_idx = perm[1]
        outboard_idx = perm[2]

        df = (;
            time_s=results.time,
            tip_displacement_m=results.tip_displacement,
            tip_twist_deg=rad2deg.(results.tip_twist),
            inboard_pitch_deg=rad2deg.(results.propeller_pitch[inboard_idx]),
            inboard_yaw_deg=rad2deg.(results.propeller_yaw[inboard_idx]),
            inboard_attachment_node=fill(attachment_nodes[inboard_idx],length(results.time)),
            inboard_attachment_wing_twist_deg=rad2deg.(results.propeller_attachment_twist[inboard_idx]),
            outboard_pitch_deg=rad2deg.(results.propeller_pitch[outboard_idx]),
            outboard_yaw_deg=rad2deg.(results.propeller_yaw[outboard_idx]),
            outboard_attachment_node=fill(attachment_nodes[outboard_idx],length(results.time)),
            outboard_attachment_wing_twist_deg=rad2deg.(results.propeller_attachment_twist[outboard_idx])
        )
        csvName = "multi_prop_A_C.csv"
        txtName = "multi_prop_A_C.txt"
        header = "Chang inboard/outboard response history"
        printlnText = "time_s\tinboard_pitch_deg\tinboard_yaw_deg\tinboard_attachment_node\tinboard_attachment_wing_twist_deg\toutboard_pitch_deg\toutboard_yaw_deg\toutboard_attachment_node\toutboard_attachment_wing_twist_deg\ttip_displacement_m\ttip_twist_deg"
    else
        error("This export helper supports 1 or 2 propellers. Got $(length(eta)).")
    end

    outputDir = joinpath(config.output.directory,"multi_prop")
    mkpath(outputDir)

    csv_path = joinpath(outputDir,csvName)
    open(csv_path,"w") do io
        println(io,join(keys(df),","))
        for i in eachindex(df.time_s)
            println(io,join(map(column -> column[i],values(df)),","))
        end
    end

    txt_path = joinpath(outputDir,txtName)
    open(txt_path,"w") do io
        println(io,header)
        if length(eta) == 2
            println(io,"inboard_propeller_index = $inboard_idx")
            println(io,"outboard_propeller_index = $outboard_idx")
            println(io,"inboard_attachment_node = $(attachment_nodes[inboard_idx])")
            println(io,"outboard_attachment_node = $(attachment_nodes[outboard_idx])")
            println(io,"eta_inboard = $(eta[inboard_idx])")
            println(io,"eta_outboard = $(eta[outboard_idx])")
        else
            println(io,"propeller_index = 1")
            println(io,"attachment_node = $(attachment_nodes[1])")
            println(io,"eta = $(eta[1])")
        end
        println(io,"")
        println(io,printlnText)
        for i in eachindex(df.time_s)
            if length(eta) == 2
                println(
                    io,
                    join(
                        [
                            df.time_s[i],
                            df.inboard_pitch_deg[i],
                            df.inboard_yaw_deg[i],
                            df.inboard_attachment_node[i],
                            df.inboard_attachment_wing_twist_deg[i],
                            df.outboard_pitch_deg[i],
                            df.outboard_yaw_deg[i],
                            df.outboard_attachment_node[i],
                            df.outboard_attachment_wing_twist_deg[i],
                            df.tip_displacement_m[i],
                            df.tip_twist_deg[i]
                        ],"\t")
                )
            else
                println(
                    io,
                    join(
                        [
                            df.time_s[i],
                            df.propeller_1_pitch_deg[i],
                            df.propeller_1_yaw_deg[i],
                            df.propeller_1_attachment_node[i],
                            df.propeller_1_attachment_wing_twist_deg[i],
                            df.tip_displacement_m[i],
                            df.tip_twist_deg[i]
                        ],"\t")
                )
            end
        end
    end

    println("Organized history CSV written to: $csv_path")
    println("Organized history text file written to: $txt_path")
    return (; csv_path,txt_path,df)
end

# Inspect chang_run.model, chang_run.workspace, or chang_run.solution as needed.
#if abspath(PROGRAM_FILE) == @__FILE__
chang_run = run_chang_example()
results = chang_run.results
organized = export_chang_organized_history(chang_run)
#end
