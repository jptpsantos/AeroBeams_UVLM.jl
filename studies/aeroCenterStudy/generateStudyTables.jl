using AeroBeams

# Tables derived from the static aeroelastic VLM solutions of the Pazy wing (test/referenceData/PazyVLM) for the study of the aerodynamic center and of the tangent lift slope (see README.md). Same layout, anchors and rules as dev/Pazy/General/generateSweptPazyTipLossFunctionDataDeformed.jl (one row per sweep angle, root angle of attack and airspeed, followed by the values at the span stations):
# - aeroCenterShiftVLM.csv: shift of the aerodynamic center of the sections relative to the quarter-chord, Δx = x_cp - 1/4 = -cm_c4/cn (chord-normalized, nose-up positive cm about the quarter-chord), at the static equilibrium of the VLM
# - tipLossTangentVLM.csv: tangent lift slope of the sections, normalized by that of the airfoil (the tip loss function of the study of the tangent slope), from the fixed-shape perturbations of the VLM (Δcn per change of local incidence)
vlmDir = joinpath(pkgdir(AeroBeams), "test", "referenceData", "PazyVLM")
outDir = @__DIR__

# Speed of sound used for the compressibility factor [m/s] and compressibility factor
a = standard_atmosphere(0).a
β_of(U) = sqrt(1-(U/a)^2)

# Limits of the shift of the aerodynamic center, as cm/cn is ill-conditioned where the load of the section vanishes (the tip of the bent wings). Limits of the tangent lift slope that are considered valid (the fixed-shape perturbations are ill-conditioned where the local incidence does not change)
Δxmin, Δxmax = -0.15, 0.05
Wmin, Wmax = 0.05, 1.3

# Read the data of all conditions (sweep angle [deg], root angle of attack [deg], airspeed [m/s])
Δx = Dict{NTuple{3,Float64},Vector{Float64}}()
Wtan = Dict{NTuple{3,Float64},Vector{Float64}}()
Wtanvalid = Dict{NTuple{3,Float64},Vector{Bool}}()
Wsec = Dict{NTuple{3,Float64},Vector{Float64}}()
invalid = NTuple{3,Float64}[]
for aoaDir in filter(d -> occursin(r"^aoa[\dp]+deg$", d), readdir(vlmDir))
    readcsv(file) = begin
        data,header = readdlm(joinpath(vlmDir, aoaDir, file), ',', header=true)
        return name -> data[:, findfirst(==(name), vec(header))]
    end
    cnData, cmData, pertData = readcsv("spanwise_cn_alpha.csv"), readcsv("spanwise_cm.csv"), readcsv("spanwise_cn_alpha_perturbations.csv")
    θ = only(unique(round.(cnData("alpha_deg"), digits=6)))
    for case in unique(cnData("case"))
        rows = cnData("case") .== case
        Λ = only(unique(round.(cnData("sweep_deg")[rows], digits=6)))
        U = only(unique(round.(cnData("U_m_s")[rows], digits=6)))
        # The three files have the same span stations of each case
        rowsCm, rowsPert = cmData("case") .== case, (pertData("case") .== case) .& (pertData("mode") .== "fixed_shape")
        @assert cmData("eta")[rowsCm] ≈ cnData("eta")[rows] ≈ pertData("eta")[rowsPert] "the span stations of the files of $aoaDir (case $case) differ"
        if !all(cnData("central_converged")[rows] .== true) || !all(cnData("cn_alpha_local_valid")[rows] .== true)
            push!(invalid, (Λ,θ,U))
            continue
        end
        cn = Float64.(cnData("cn")[rows])
        cmc4 = Float64.(cmData("cm_quarter_chord")[rowsCm])
        Δx[(Λ,θ,U)] = clamp.(-cmc4 ./ cn, Δxmin, Δxmax)
        Wsec[(Λ,θ,U)] = Float64.(cnData("cn_alpha_local_per_rad")[rows]) * β_of(U)/(2π)
        Wtan[(Λ,θ,U)] = Float64.(pertData("cn_alpha_local_perturbation_per_rad")[rowsPert]) * β_of(U)/(2π)
        Wtanvalid[(Λ,θ,U)] = pertData("valid")[rowsPert] .== true
    end
end
ΛRange = sort(unique(vcat(getindex.(keys(Δx), 1), getindex.(invalid, 1))))
θData = sort(unique(vcat(getindex.(keys(Δx), 2), getindex.(invalid, 2))))
URange = sort(unique(vcat(getindex.(keys(Δx), 3), getindex.(invalid, 3))))
θAnchor, UAnchor = first(θData), first(URange)
@assert θAnchor <= 1e-2 && UAnchor <= 0.1 "the VLM data must have a condition with root angle of attack and airspeed close to zero, as the anchor of the tables"
θRange = vcat(0.0, θData[2:end])

# Conditions in which the VLM does not converge: use the data of the closest airspeed with data at the same sweep and root angle
for (Λ,θ,U) in invalid
    UValid = [k[3] for k in keys(Δx) if k[1] == Λ && k[2] == θ]
    UClosest = UValid[argmin(abs.(UValid .- U))]
    for d in (Δx, Wtan, Wsec, Wtanvalid)
        d[(Λ,θ,U)] = d[(Λ,θ,UClosest)]
    end
    println("Using the data at U=$UClosest m/s for Λ=$(Λ)°, θ=$(θ)°, U=$U m/s")
end

# The undeformed wing (zero root angle) has the values of the anchor at all airspeeds (the tangent lift slope of the undeformed wing is the secant one, with the compressibility factor of each airspeed)
for Λ in ΛRange, U in URange
    Δx[(Λ,0.0,U)] = Δx[(Λ,θAnchor,UAnchor)]
    Wsec[(Λ,0.0,U)] = Wsec[(Λ,θAnchor,UAnchor)] * β_of(U)/β_of(UAnchor)
    Wtan[(Λ,0.0,U)] = Wsec[(Λ,0.0,U)]
    Wtanvalid[(Λ,0.0,U)] = trues(length(Wtan[(Λ,0.0,U)]))
end

# Tangent lift slope: where it is not valid (or out of the limits), use that of the undeformed wing
nReplaced = 0
for Λ in ΛRange, θ in θRange, U in URange
    bad = .!Wtanvalid[(Λ,θ,U)] .| .!isfinite.(Wtan[(Λ,θ,U)]) .| (Wtan[(Λ,θ,U)] .< Wmin) .| (Wtan[(Λ,θ,U)] .> Wmax)
    Wtan[(Λ,θ,U)] = ifelse.(bad, Wsec[(Λ,0.0,U)], Wtan[(Λ,θ,U)])
    global nReplaced += count(bad)
end
println("Entries of the tangent lift slope replaced by that of the undeformed wing: $nReplaced")

# Save the tables
nStations = length(first(values(Δx)))
header = ["sweep_deg" "aoa_deg" "U_m_s" ["s=$(k-1)/$(nStations-1)" for k in 1:nStations]...]
for (file, conditions) in (("aeroCenterShiftVLM.csv", Δx), ("tipLossTangentVLM.csv", Wtan))
    table = [[Λ θ U transpose(conditions[(Λ,θ,U)])] for Λ in ΛRange for θ in θRange for U in URange]
    open(joinpath(outDir, file), "w") do io
        writedlm(io, header, ',')
        writedlm(io, vcat(table...), ',')
    end
end
println("Finished generateStudyTables.jl")
