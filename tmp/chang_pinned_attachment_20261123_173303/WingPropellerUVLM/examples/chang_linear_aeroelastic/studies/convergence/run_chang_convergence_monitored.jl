# Run from PowerShell with --threads=10 and the local studies project.
using Dates,DelimitedFiles,Printf

const ENTRY = joinpath(@__DIR__,"chang_convergence_aeroelastic.jl")
const OUTPUT = normpath(joinpath(@__DIR__,"..","..","output"))
const RATE_COLUMNS = (
    "pitch_lambda_per_s" => "pitch lambda",
    "yaw_lambda_per_s" => "yaw lambda",
    "wing_bending_lambda_per_s" => "bending lambda",
    "wing_torsion_lambda_per_s" => "torsion lambda"
)

function read_campaign_table(path)
    rows,header = readdlm(path,','; header=true)
    names = string.(vec(header))
    return [Dict(name => rows[i,j] for (j,name) in enumerate(names)) for i in axes(rows,1)]
end

# The solver rewrites its tables after each case.
function campaign_snapshot(path)
    isfile(path) || return Dict{String,Any}[]
    try
        return read_campaign_table(path)
    catch exception
        exception isa InterruptException && rethrow()
        return Dict{String,Any}[]
    end
end

reused(row) = string(row["reused"]) == "true"

function campaign_log(io,message)
    line = "$(Dates.format(now(),"yyyy-mm-dd HH:MM:SS")) $message"
    println(io,line)
    flush(io)
    println(line)
    flush(stdout)
    return nothing
end

function live_case_stage(directory,case)
    path = joinpath(directory,case["label"],"run.log")
    isfile(path) || return "starting | $(case["steps"]) planned steps"
    tail = open(path) do io
        seekend(io)
        seek(io,max(0,position(io)-4096))
        read(io,String)
    end
    if occursin("Simulation finished",tail) || occursin("VALIDATION_PASSED",tail)
        return "solver finished | processing damping fits"
    end
    steps = collect(eachmatch(r"(?m)^Step (\d+)/(\d+)",tail))
    isempty(steps) && return "initializing | $(case["steps"]) planned steps"
    done,total = parse.(Int,last(steps).captures)
    return @sprintf("step %d/%d (%.1f%%)",done,total,100done/total)
end

function monitor_campaign(io,directory,worker;dryRun=false)
    planPath = joinpath(directory,"runtime_plan.csv")
    matrixPath = joinpath(directory,"case_matrix.csv")
    resultPath = joinpath(directory,"convergence.csv")
    plan = campaign_snapshot(planPath)
    matrix = Dict(r["label"] => r for r in campaign_snapshot(matrixPath))
    while process_running(worker) && (isempty(plan) || any(r -> !haskey(matrix,r["label"]),plan))
        sleep(1)
        plan = campaign_snapshot(planPath)
        matrix = Dict(r["label"] => r for r in campaign_snapshot(matrixPath))
    end
    isempty(plan) && return
    n = length(plan)
    uniqueTotal = count(r -> !reused(r),plan)
    labels = Dict(r["label"] => i for (i,r) in enumerate(plan))
    campaign_log(io,"PLAN | $n entries | $uniqueTotal unique solves")
    for (i,row) in enumerate(plan)
        label,family = row["label"],row["family"]
        setting = get(get(matrix,label,Dict()),family,"")
        reuse = reused(row) ? " | reused solve" : ""
        campaign_log(io,"PLAN $i/$n $label | $family=$setting | $(row["steps"]) steps | $(row["wake_rows"]) wake rows$reuse")
    end
    dryRun && return
    outcomes = Dict{String,Dict}()
    lastProgress = ""
    while true
        running = process_running(worker)
        for row in campaign_snapshot(resultPath)
            label = string(get(row,"label",""))
            haskey(labels,label) && !haskey(outcomes,label) && get(row,"elapsed_s",nothing) isa Real ||
                continue
            outcomes[label] = row
            rates = [@sprintf("%s=%.6g s^-1",name,row[key]) for (key,name) in RATE_COLUMNS
                if get(row,key,nothing) isa Real]
            detail = row["status"] == "failed" ? string(row["reason"]) : join(rates," | ")
            reuse = reused(row) ? " | reused" : ""
            campaign_log(io,"RESULT $(labels[label])/$n $label | $(row["status"]) | $(@sprintf("%.1f",row["elapsed_s"])) s$reuse | $detail")
        end
        uniqueDone = count(r -> !reused(r),values(outcomes))
        counts = "$(length(outcomes))/$n entries | $uniqueDone/$uniqueTotal unique solves"
        if !running
            statuses = sort(unique(r["status"] for r in values(outcomes)))
            result = join(["$s=$(count(r -> r["status"] == s,values(outcomes)))" for s in statuses],", ")
            campaign_log(io,"EXIT | $counts | $result")
            return nothing
        end
        next = findfirst(r -> !reused(r) && !haskey(outcomes,r["label"]),plan)
        progress = if isnothing(next)
            "PROGRESS $counts | finalizing the campaign"
        else
            row = plan[next]
            "PROGRESS $counts | [$next/$n] $(row["label"]) | $(live_case_stage(directory,row))"
        end
        if progress != lastProgress
            campaign_log(io,progress)
            lastProgress = progress
        end
        sleep(5)
    end
end

function run_monitored(;dryRun=false,dampingTarget=:propeller)
    dampingTarget in (:wing,:propeller) || throw(ArgumentError("Use --damping-target=wing or --damping-target=propeller"))
    mkpath(OUTPUT)
    prefix = "convergence_aeroelastic_$(Threads.nthreads())threads_$(Dates.format(now(),"yyyymmdd_HHMMSS"))_"
    directory = mktempdir(OUTPUT;prefix,cleanup=false)
    out = joinpath(directory,"campaign.stdout.log")
    err = joinpath(directory,"campaign.stderr.log")
    live = joinpath(directory,"campaign.live.log")
    project = dirname(Base.active_project())
    args = ["--threads=$(Threads.nthreads())","--project=$project",@__FILE__,"--worker",directory,
        string(dryRun),string(dampingTarget)]
    cmd = Cmd([Base.julia_cmd().exec; args])
    open(live,"w") do log
        campaign_log(log,"Chang convergence | $(Threads.nthreads()) Julia threads | damping target $dampingTarget | $directory")
        worker = run(pipeline(cmd;stdout=out,stderr=err);wait=false)
        campaign_log(log,"Julia PID $(getpid(worker)) | initializing study")
        monitor_campaign(log,directory,worker;dryRun)
        wait(worker)
        if !success(worker)
            campaign_log(log,"Solver exited with code $(worker.exitcode). See $err")
            exit(worker.exitcode)
        end
    end
    return directory
end

function monitored_cli(args=ARGS)
    dryRun = false
    dampingTarget = :propeller
    for arg in args
        if arg == "--dry-run"
            dryRun = true
        elseif startswith(arg,"--damping-target=")
            dampingTarget = Symbol(last(split(arg,'=';limit=2)))
        else
            throw(ArgumentError("Unknown argument $arg. Use --dry-run or --damping-target=wing|propeller."))
        end
    end
    return run_monitored(;dryRun,dampingTarget)
end

if abspath(PROGRAM_FILE) == @__FILE__
    if !isempty(ARGS) && first(ARGS) == "--worker"
        get!(ENV,"GKSwstype","100")
        include(ENTRY)
        const C = ChangConvergenceAeroelastic
        dampingTarget = length(ARGS) >= 4 ? Symbol(ARGS[4]) : :propeller
        settings = C.create_study(; runMode=:full,dryRun=parse(Bool,ARGS[3]),dampingTarget).settings
        study = C.create_study(settings; outputDirectory=ARGS[2])
        C.WingPropellerUVLMStudies.solve!(study)
    else
        monitored_cli()
    end
end
