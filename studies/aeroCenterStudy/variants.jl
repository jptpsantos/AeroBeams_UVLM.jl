# Variants of the aerodynamic model of the Pazy wing for the flutter and static scripts (see README.md and runVariant.jl). Included by runVariant.jl after it defines `key` (the script) and `OVERRIDES` (text replacements of the script). The variants are selected with environment variables:
#   AF       airfoil: script (as in the script) | flat | pazy (NACA0018Pazy) | xac:<x> (constant x_ac, cnα=2π) | xacU:<U1>,<x1>;<U2>,<x2>;... (x_ac as a function of airspeed [m/s], cnα=2π) | timmerx (x_ac of Timmer's data, cnα=2π) | timmerfull (NACA0018Timmer) | timmerfull:<Δ> (NACA0018Timmer with Δ added to x_ac)
#   FWL      1/0: finite-wing lag on/off (the default of Indicial() is on)
#   WS       constant factor of the VLM tip loss function
#   TLMODE   pkg (default: the :VLM_tan_fixed tip loss of the package, tangent lift slope since 2026-11-21) | sec (secant lift slope, tipLossSecantVLM.csv, the :VLM table before that date) | tan (tangent lift slope, tipLossTangentVLM.csv, same as the package) | swap (secant for the steady solution, tangent for the eigenproblem)
#   XS       1: shift of the aerodynamic center along the span from aeroCenterShiftVLM.csv, XSSCALE multiplies it (requires the prototype of aeroCenterShift.patch)
#   WN       1: tip loss function scaled by the ratio of the lift slope of the Nastran model of Technion to that of the VLM (θ=5°; Λ=0 and 20°; U ≤ 60 m/s, held at higher airspeeds)
#   TLTH0    1: tip loss function of the undeformed wing (θ=0) at any root pitch angle
using DelimitedFiles, LinearAlgebra, StaticArrays, LinearInterpolations

# Airfoil with the given Mach numbers (default: constant) and aerodynamic center locations, and the lift slope of the flat plate
function airfoil_with_xac(MaRng, xac; cnα=2π, cd₀=1e-2, name="X")
    attached = isnothing(MaRng) ? (cd₀=cd₀, cm₀=0.0, cnα=cnα, x_ac=xac) : (MaRng=MaRng, cd₀=cd₀, cm₀=0.0, cnα=cnα, x_ac=xac)
    bl = isnothing(MaRng) ? (cd₀=cd₀, cm₀=0.0, cnα=cnα, K₀=0.25 .- xac) : (MaRng=MaRng, cd₀=cd₀, cm₀=0.0, cnα=cnα, K₀=0.25 .- xac)
    return create_Airfoil(database=modified_database("NACA0018"; attached=attached, bli=bl, blo=bl, newName=name))
end

# Mach numbers and aerodynamic center locations of the data of Timmer (2008), as in NACA0018Timmer (chord of the straight Pazy wing, sea level)
function timmer_xac()
    atm = standard_atmosphere(0)
    ReRng = [0.15, 0.30, 0.50, 0.70, 1.00]*1e6
    MaRng = ReRng * atm.μ / (atm.ρ * atm.a * 0.0988504)
    cl₇ = [0.830, 0.761, 0.693, 0.691, 0.681]
    cm₇ = [-0.0064, 0.0060, 0.0197, 0.0215, 0.0247]
    return MaRng, 0.25 .- cm₇./cl₇
end

function make_airfoil(spec)
    atm = standard_atmosphere(0)
    if spec == "flat"
        return deepcopy(FlatPlate)
    elseif spec == "pazy"
        return deepcopy(NACA0018Pazy)
    elseif startswith(spec, "xac:")
        return airfoil_with_xac(nothing, parse(Float64, spec[5:end]))
    elseif startswith(spec, "xacU:")
        pairs_ = [parse.(Float64, split(p, ",")) for p in split(spec[6:end], ";")]
        return airfoil_with_xac([p[1] for p in pairs_]/atm.a, [p[2] for p in pairs_])
    elseif spec == "timmerx"
        MaRng, xac = timmer_xac()
        return airfoil_with_xac(MaRng, xac)
    elseif spec == "timmerfull"
        return deepcopy(NACA0018Timmer)
    elseif startswith(spec, "timmerfull:")
        MaRng, xac = timmer_xac()
        return create_Airfoil(database=modified_database("NACA0018Timmer"; attached=(MaRng=MaRng, x_ac=xac .+ parse(Float64, spec[12:end])), newName="NACA0018TimmerOffset"))
    end
    error("unknown airfoil variant $spec")
end

AF = get(ENV, "AF", "script")
if AF != "script"
    push!(OVERRIDES, (key == "s20" ? "airfoil = deepcopy(FlatPlate)" : "airfoil = deepcopy(NACA0018Pazy)") => "airfoil = make_airfoil(\"$AF\")")
end
if haskey(ENV, "FWL")
    push!(OVERRIDES, "aeroSolver = Indicial()" => "aeroSolver = Indicial(finiteWingLag=$(ENV["FWL"] == "1"))")
end
if get(ENV, "TLTH0", "0") == "1"
    push!(OVERRIDES, "beamWingtipRod=beamWingtipRod," => "beamWingtipRod=beamWingtipRod,tipLossPitchAngle=0.0,")
end

## Tip loss function: a single replacement of AeroBeams.tip_loss_vlm
const TL_MODE = Ref(:sec)
function load_table(file)
    data = readdlm(file, ',', Float64; skipstart=1)
    Λ, θ, U = sort(unique(data[:,1])), sort(unique(data[:,2])), sort(unique(data[:,3]))
    nS = size(data,2) - 3
    D = Array{SVector{nS,Float64}}(undef, length(Λ), length(θ), length(U))
    for row in eachrow(data)
        D[searchsortedfirst(Λ,row[1]), searchsortedfirst(θ,row[2]), searchsortedfirst(U,row[3])] = SVector{nS,Float64}(row[4:end])
    end
    return LinearInterpolations.Interpolate((deg2rad.(Λ), deg2rad.(θ), U), D; extrapolate=:replicate)
end
# Ratio of the lift slope distribution of the Nastran model of Technion (θ=5°) to that of the VLM, at the span stations of the Nastran data, as a function of airspeed (10 to 60 m/s) and sweep angle (0 and 20°, linear in between)
function nastran_ratio()
    refdir = joinpath(pkgdir(AeroBeams), "test", "referenceData", "sweptPazy")
    ratio = Dict{Int,Any}()
    for (Λd,file) in ((0,"clalapha_Lambda0_aoa5_NastranTechnion.txt"), (20,"clalapha_Lambda20_aoa5_NastranTechnion.txt"))
        ref = readdlm(joinpath(refdir, file))
        rows = Vector{Vector{Float64}}(); sN0 = Float64[]
        for j in 1:6
            sN, cN = ref[2j-1,:], ref[2j,:]
            ok = .!isnan.(sN) .& .!isnan.(cN)
            sN, cN = sN[ok], cN[ok]
            p = sortperm(sN); sN, cN = sN[p], cN[p]
            W = AeroBeams.tip_loss_function_Pazy(:VLM_tan_fixed, deg2rad(5), 10.0*j, deg2rad(Λd))
            push!(rows, cN ./ (2π .* W.(sN))); sN0 = sN
        end
        ratio[Λd] = (sN0, rows)
    end
    return ratio
end
function nastran_ratio_fun(U, Λ)
    Uc = clamp(U, 10.0, 60.0); j = clamp(floor(Int, Uc/10), 1, 5); w = (Uc - 10j)/10
    function atΛ(Λd)
        sN, rows = NASTRAN_RATIO[Λd]
        r = (1-w) .* rows[j] .+ w .* rows[j+1]
        return s -> s <= sN[1] ? r[1] : s >= sN[end] ? r[end] : begin
            k = searchsortedlast(sN, s); r[k] + (r[k+1]-r[k])*(s-sN[k])/(sN[k+1]-sN[k])
        end
    end
    a = clamp(rad2deg(Λ)/20, 0.0, 1.0); f0, f20 = atΛ(0), atΛ(20)
    return s -> (1-a)*f0(s) + a*f20(s)
end
useNastran = get(ENV, "WN", "0") == "1"
useNastran && (global NASTRAN_RATIO = nastran_ratio())
WSCALE = parse(Float64, get(ENV, "WS", "1.0"))
tlmode = get(ENV, "TLMODE", "pkg")
useShift = get(ENV, "XS", "0") == "1"
xsscale = parse(Float64, get(ENV, "XSSCALE", "1.0"))
if useShift
    @assert isdefined(AeroBeams, :TipLossWithShift) "the shift of the aerodynamic center requires the prototype of aeroCenterShift.patch: run with --project=<copy of the package with the patch applied>"
    global XS_INTERP = load_table(joinpath(@__DIR__, "aeroCenterShiftVLM.csv"))
end
tlmode in ("tan", "swap") && (global XT_INTERP = load_table(joinpath(@__DIR__, "tipLossTangentVLM.csv")))
tlmode in ("sec", "swap") && (global XS_SEC = load_table(joinpath(@__DIR__, "tipLossSecantVLM.csv")))
if WSCALE != 1.0 || useShift || useNastran || tlmode != "pkg"
    Core.eval(AeroBeams, quote
        function tip_loss_vlm(type,θ,U,Λ)
            Ws = $(tlmode in ("sec", "swap")) ? Main.XS_SEC((Λ,θ,U)) : _VLM_tipLoss[type]((Λ,θ,U))
            Wt = $(tlmode in ("tan", "swap")) ? Main.XT_INTERP((Λ,θ,U)) : Ws
            R = $(useNastran) ? Main.nastran_ratio_fun(U,Λ) : nothing
            wfun = s -> begin
                W = $(tlmode == "tan") ? Wt : ($(tlmode == "swap") ? (Main.TL_MODE[] == :tan ? Wt : Ws) : Ws)
                n = length(W); x = clamp(s,0,1)*(n-1); i = min(floor(Int,x)+1,n-1)
                v = $(WSCALE) * (@inbounds W[i] + (x-(i-1))*(W[i+1]-W[i]))
                return $(useNastran) ? v*R(s) : v
            end
            $(useShift) || return wfun
            X = $(useShift) ? Main.XS_INTERP((Λ,θ,U)) : nothing
            dfun = s -> begin
                n = length(X); x = clamp(s,0,1)*(n-1); i = min(floor(Int,x)+1,n-1)
                @inbounds $(xsscale)*(X[i] + (x-(i-1))*(X[i+1]-X[i]))
            end
            return TipLossWithShift(wfun, dfun)
        end
    end)
end

# Eigenproblem with the tangent lift slope around the steady solution of the secant one: the steady loads are those of the secant lift slope (static equilibrium of the VLM), while the perturbations of the eigenproblem (flutter) see the tangent lift slope of the VLM
function solve_swap!(problem)
    AeroBeams.precompute_motion_basis_A!(problem); AeroBeams.precompute_BCs_history!(problem)
    AeroBeams.precompute_distributed_loads!(problem); AeroBeams.precompute_flap_deflection!(problem)
    AeroBeams.update_motion_basis_A!(problem,1); AeroBeams.update_BCs_state!(problem,1)
    TL_MODE[] = :sec
    AeroBeams.solve_steady!(problem)
    TL_MODE[] = :tan
    AeroBeams.foreach_element(problem.model.assemblyData.elementGroups) do element
        if !isnothing(element.aero.surface)
            ld = element.aero.loadsDistribution
            ld.ϖMid = ld.ϖ(ld.Uref, 0.5)
        end
    end
    AeroBeams.j!(problem)
    problem.systemArrays.luCache = lu(problem.systemArrays.jacobian; check=false)
    AeroBeams.i!(problem)
    AeroBeams.solve_eigen!(problem)
    TL_MODE[] = :sec
    return nothing
end
tlmode == "swap" && push!(OVERRIDES, "solve!(problem)" => "solve_swap!(problem)")
