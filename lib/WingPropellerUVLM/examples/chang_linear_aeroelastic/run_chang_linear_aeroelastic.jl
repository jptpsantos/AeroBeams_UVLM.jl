# Chang wing–propeller aeroelastic response.
# Edit chang_case.jl for the case, numerical controls, and output.
# CHANG_* environment variables override those defaults when config is loaded.

import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")))

if !isdefined(@__MODULE__, :ChangAeroelastic) ||
        !isdefined(ChangAeroelastic, :load_chang_configuration)
    include(joinpath(@__DIR__, "src", "ChangAeroelastic.jl"))
end
using .ChangAeroelastic
using CSV
using DataFrames

function run_chang_example()
    config = load_chang_configuration()
    return run_chang(config)
end

function export_chang_organized_history(chang_run)
    config = chang_run.config
    results = chang_run.results
    eta = config.propeller.attachment_eta

    if length(eta) == 1
        df = DataFrame(
            time_s = results.time,
            tip_displacement_m = results.tip_displacement,
            tip_twist_deg = rad2deg.(results.tip_twist),
            propeller_1_pitch_deg = rad2deg.(results.propeller_pitch[1]),
            propeller_1_yaw_deg = rad2deg.(results.propeller_yaw[1]),
        )
        csv_name = "single_propeller_history.csv"
        txt_name = "single_propeller_history.txt"
        header = "Chang single-propeller response history"
        println_text = "time_s\tpropeller_1_pitch_deg\tpropeller_1_yaw_deg\ttip_displacement_m\ttip_twist_deg"
    elseif length(eta) == 2
        perm = sortperm(eta)
        inboard_idx = perm[1]
        outboard_idx = perm[2]

        df = DataFrame(
            time_s = results.time,
            tip_displacement_m = results.tip_displacement,
            tip_twist_deg = rad2deg.(results.tip_twist),
            inboard_pitch_deg = rad2deg.(results.propeller_pitch[inboard_idx]),
            inboard_yaw_deg = rad2deg.(results.propeller_yaw[inboard_idx]),
            outboard_pitch_deg = rad2deg.(results.propeller_pitch[outboard_idx]),
            outboard_yaw_deg = rad2deg.(results.propeller_yaw[outboard_idx]),
        )
        csv_name = "multi_prop_A_C.csv"
        txt_name = "multi_prop_A_C.txt"
        header = "Chang inboard/outboard response history"
        println_text = "time_s\tinboard_pitch_deg\tinboard_yaw_deg\toutboard_pitch_deg\toutboard_yaw_deg\ttip_displacement_m\ttip_twist_deg"
    else
        error("This export helper supports 1 or 2 propellers. Got $(length(eta)).")
    end

    output_dir = joinpath(@__DIR__, "output", "multi_prop")
    mkpath(output_dir)

    csv_path = joinpath(output_dir, csv_name)
    CSV.write(csv_path, df)

    txt_path = joinpath(output_dir, txt_name)
    open(txt_path, "w") do io
        println(io, header)
        if length(eta) == 2
            println(io, "inboard_propeller_index = $(sortperm(eta)[1])")
            println(io, "outboard_propeller_index = $(sortperm(eta)[2])")
            println(io, "eta_inboard = $(eta[sortperm(eta)[1]])")
            println(io, "eta_outboard = $(eta[sortperm(eta)[2]])")
        else
            println(io, "propeller_index = 1")
            println(io, "eta = $(eta[1])")
        end
        println(io, "")
        println(io, println_text)
        for i in eachindex(df.time_s)
            if length(eta) == 2
                println(
                    io,
                    join([
                        df.time_s[i],
                        df.inboard_pitch_deg[i],
                        df.inboard_yaw_deg[i],
                        df.outboard_pitch_deg[i],
                        df.outboard_yaw_deg[i],
                        df.tip_displacement_m[i],
                        df.tip_twist_deg[i],
                    ], "\t"),
                )
            else
                println(
                    io,
                    join([
                        df.time_s[i],
                        df.propeller_1_pitch_deg[i],
                        df.propeller_1_yaw_deg[i],
                        df.tip_displacement_m[i],
                        df.tip_twist_deg[i],
                    ], "\t"),
                )
            end
        end
    end

    println("Organized history CSV written to: $csv_path")
    println("Organized history text file written to: $txt_path")
    return (; csv_path, txt_path, df)
end

# Inspect chang_run.model, chang_run.workspace, or chang_run.solution as needed.
#if abspath(PROGRAM_FILE) == @__FILE__
    chang_run = run_chang_example()
    results = chang_run.results
    organized = export_chang_organized_history(chang_run)
#end
