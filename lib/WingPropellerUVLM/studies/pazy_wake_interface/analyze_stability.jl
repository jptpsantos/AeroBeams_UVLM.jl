# Run with the repository's examples project in a separate Julia process.
using AeroBeams, DelimitedFiles, LinearAlgebra, Printf

const example = normpath(joinpath(@__DIR__,"..","..","examples","pazy_aerobeams_coupled"))
const output = mkpath(joinpath(@__DIR__,"output"))
const pitch = deg2rad(1.0)
const profile_type = "FixedExponential-99999"
const profile = Ref{Tuple{Vector{Float64},Vector{Float64}}}()
BLAS.set_num_threads(1)

function interpolate_weight(eta,values,s)
    i = clamp(searchsortedlast(eta,s),1,length(eta)-1)
    t = (s - eta[i]) / (eta[i+1] - eta[i])
    return (1-t) * values[i] + t * values[i+1]
end

# Specialize only the standalone study's tip-loss selection. create_Pazy
# constructs the aerodynamic elements with this function from the start.
function AeroBeams.tip_loss_function_Pazy(type::String,θ::Real,U::Float64,Λ::Real)
    type == profile_type || return invoke(AeroBeams.tip_loss_function_Pazy,
        Tuple{String,Real,Real,Real},type,θ,U,Λ)
    eta,values = profile[]
    return s -> interpolate_weight(eta,values,s)
end

function read_profile(path)
    data,header = readdlm(path,',',Float64; header=true)
    col(name) = only(findall(==(name),vec(header)))
    rows = data[:,col("U_m_s")] .== 10.0
    eta = data[rows,col("eta")]
    values = data[rows,col("cn_alpha_fixed_shape_local_per_rad")] ./ (2π)
    @assert length(eta) == 21 && issorted(eta) && first(eta) == 0.0 && last(eta) == 1.0
    @assert all(isfinite,values)
    return eta,values
end

function solve_modal_sweep(speeds)
    eigensolutions = EigenSolution[]
    tip_deflections = Float64[]
    x0 = Float64[]
    _,span,chord,spar,_ = geometrical_properties_Pazy(; mesh=:refined)
    for U in speeds
        model,_ = create_Pazy(aeroSolver=Indicial(),airfoil=deepcopy(FlatPlate),mesh=:refined,
            hasTipCorrection=true,tipLossType=profile_type,upright=true,pitchAngle=0.0,
            flightDirection=[sin(pitch),cos(pitch),0.0],gravityVector=zeros(3),airspeed=U)
        problem = create_EigenProblem(; model,eigenmethod=:krylovkit,nev=65,x0)
        solve!(problem)
        problem.systemSolver.loadStepState.convergedFinalSolution || error("Trim failed at U=$U")
        push!(eigensolutions,problem.eigensolution)
        tip = last(problem.loadstepHistory.states)[length(model.beams[1].elements)].nodal
        twist = twist_angle(tip.p2_b; outputUnit=:deg)
        push!(tip_deflections,-(tip.u2[1] - chord*(0.5-spar)*sind(twist))/span)
        x0 = copy(problem.systemArrays.x)
        isinteger(U) && iszero(U % 10) && @printf("Solved U=%g m/s\n",U)
        flush(stdout)
    end
    frequencies,growth,_ = mode_tracking(eigensolutions)
    rads2Hz!(frequencies)
    frequencies,growth,_ = filter_tracked_modes(frequencies,growth; fmin=0.1,maxTracks=5)
    return frequencies,growth,tip_deflections
end

function export_stability(tag,speeds,frequencies,growth,tip_deflections)
    results = find_flutter(speeds,frequencies,growth,tip_deflections)
    onsets,onsetFreqs,onsetModes,_,offsets,offsetFreqs,offsetModes,_ = results
    header = ["kind" "mode" "speed_m_s" "frequency_Hz"]
    onrows = hcat(fill("onset",length(onsets)),onsetModes,onsets,onsetFreqs)
    offrows = hcat(fill("offset",length(offsets)),offsetModes,offsets,offsetFreqs)
    writedlm(joinpath(output,"stability_$(tag)_boundaries.csv"),vcat(header,onrows,offrows),',')
    header = hcat(["U_m_s" "tip_deflection_over_span"],
        permutedims(["frequency_mode_$i" for i in axes(frequencies,1)]),
        permutedims(["growth_mode_$i" for i in axes(growth,1)]))
    writedlm(joinpath(output,"stability_$(tag)_modes.csv"),
        vcat(header,hcat(speeds,tip_deflections,transpose(frequencies),transpose(growth))),',')
    println("$tag crossings: ",vcat(onrows,offrows))
end

function run_stability(tag)
    sources = Dict(
        "baseline" => joinpath(example,"output","static_spanwise_cn_cm","aoa1deg","spanwise_cn_alpha_clean.csv"),
        "corrected" => joinpath(output,"coupled","aoa1deg","spanwise_cn_alpha_clean.csv"))
    profile[] = read_profile(sources[tag])
    speeds = vcat(0.001,10.0:10.0:40.0,45.0:1.0:54.0,55.0:0.5:85.0,
        86.0:1.0:95.0,100.0:10.0:120.0)
    frequencies,growth,tip_deflections = solve_modal_sweep(speeds)
    export_stability(tag,speeds,frequencies,growth,tip_deflections)
end

run_stability(only(ARGS))
