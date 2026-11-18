using AeroBeams

# Directory of the static aeroelastic VLM solutions of the Pazy wing: one folder per sweep angle (S<Λ>), each with one folder per root angle of attack (aoa<θ>deg, with "p" as decimal separator), which holds the spanwise distributions (spanwise_cn_alpha.csv) for a range of airspeeds
vlmDir = joinpath(pkgdir(AeroBeams), "test", "referenceData", "PazyVLM")

# Speed of sound used for the compressibility factor [m/s] (that of the airfoil's lift slope in AeroBeams, at sea level)
a = standard_atmosphere(0).a

# Compressibility factor as a function of airspeed
β_of(U) = sqrt(1-(U/a)^2)

# Read the lift slope along the span (cn_alpha_local_per_rad) and the tip out-of-plane displacement of the VLM data of all conditions (sweep angle [deg], root angle of attack [deg], airspeed [m/s])
slopes = Dict{NTuple{3,Float64},Vector{Float64}}()
tipOOP = Dict{NTuple{3,Float64},Float64}()
zeroLoad = NTuple{2,Float64}[]
eta = Float64[]
for sweepDir in filter(d -> occursin(r"^S\d+$", d), readdir(vlmDir)), aoaDir in filter(d -> occursin(r"^aoa[\dp]+deg$", d), readdir(joinpath(vlmDir, sweepDir)))
    # Read data and set the function to get its columns
    data,header = readdlm(joinpath(vlmDir, sweepDir, aoaDir, "spanwise_cn_alpha.csv"), ',', header=true)
    column(name) = data[:, findfirst(==(name), vec(header))]
    Λ = only(unique(round.(column("sweep_deg"), digits=6)))
    θ = only(unique(round.(column("alpha_deg"), digits=6)))
    # Conditions without lift (zero angle of attack) have no slope (cn/α is not defined, and the wing is undeformed): they are set further below
    valid = column("cn_alpha_local_valid") .== true
    if !any(valid)
        println("No valid lift slope (zero load) in $sweepDir/$aoaDir")
        push!(zeroLoad, (Λ,θ))
        continue
    end
    @assert all(valid) && all(column("central_converged") .== true) "VLM data of $sweepDir/$aoaDir is not valid and converged at all points"
    # Data of each airspeed
    for case in unique(column("case"))
        rows = column("case") .== case
        U = only(unique(round.(column("U_m_s")[rows], digits=6)))
        slopes[(Λ,θ,U)] = column("cn_alpha_local_per_rad")[rows]
        tipOOP[(Λ,θ,U)] = only(unique(column("tip_OOP_m")[rows]))
        # The span stations are the same for all conditions
        etaCase = column("eta")[rows]
        isempty(eta) ? append!(eta, etaCase) : @assert etaCase ≈ eta "the span stations of $sweepDir/$aoaDir differ from the others"
    end
end
URange = sort(unique(getindex.(keys(slopes), 3)))

# Without load, the wing is undeformed, and its lift slope does not depend on the airspeed (apart from the compressibility factor, below): use that of the condition closest to the undeformed state, with the lowest angle of attack and airspeed
for (Λ,θ) in zeroLoad
    θClosest = minimum(k[2] for k in keys(slopes) if k[1] == Λ)
    for U in URange
        slopes[(Λ,θ,U)] = slopes[(Λ,θClosest,URange[1])]
    end
    println("Using the lift slope at θ=$(θClosest)° and U=$(URange[1]) m/s, with tip OOP displacement of $(round(1e3*tipOOP[(Λ,θClosest,URange[1])], digits=1)) mm, for the zero load condition at Λ=$(Λ)°, θ=$(θ)°")
end

# Ranges of the conditions, which must form a complete grid
ΛRange = sort(unique(getindex.(keys(slopes), 1)))
θRange = sort(unique(getindex.(keys(slopes), 2)))
@assert length(slopes) == length(ΛRange)*length(θRange)*length(URange) "the VLM data does not cover all combinations of sweep angle, root angle of attack and airspeed"
@assert eta ≈ range(0, 1, length=length(eta)) "the span stations must be evenly spaced from the root to the tip"
println("VLM data at Λ = $ΛRange deg, θ = $θRange deg, U = $URange m/s, and $(length(eta)) span stations")

# The tip loss function is the lift slope along the span normalized by that of the airfoil, 2π/β. No sweep factor is needed: the VLM coefficient is based on the chord of the AeroBeams model and on the spanwise width projected on the unswept direction (Δℓ cosΛ), which matches the sweep correction (cosΛ) of the loads of the strips of AeroBeams. Compare the static equilibria of AeroBeams and of the VLM (tip_OOP_m and tip_twist_deg in the data) to check it
conditions = Dict(k => slope / (2π/β_of(k[3])) for (k,slope) in slopes)

# Save the tip loss function at the span stations, one row per condition
table = [[Λ θ U transpose(conditions[(Λ,θ,U)])] for Λ in ΛRange for θ in θRange for U in URange]
header = ["sweep_deg" "aoa_deg" "U_m_s" ["s=$(k-1)/$(length(eta)-1)" for k in eachindex(eta)]...]
open(joinpath(pkgdir(AeroBeams), "test", "referenceData", "sweptPazy", "tipLossFunctionVLM.csv"), "w") do io
    writedlm(io, header, ',')
    writedlm(io, vcat(table...), ',')
end

# Plot a grid condition (the data points are exactly reproduced) and a query between grid points, as used by the models. Restart Julia if the data was already loaded
Λq, θq, Uq = 15, 6, 65
Λg, θg, Ug = 20, 5, 60
s = 0:0.01:1
plt = plot(xlims=[0,1], ylims=[0,1.5], xlabel="Normalized span", ylabel="Normalized lift slope, \$c_{n_\\alpha}/(2\\pi/\\beta)\$")
scatter!(eta, conditions[(Float64(Λg),Float64(θg),Float64(Ug))], c=:black, label="VLM data: Λ=$(Λg)°, θ=$(θg)°, U=$Ug m/s")
plot!(s, tip_loss_function_Pazy("VLM-def", deg2rad(θg), Ug, deg2rad(Λg)).(s), c=:black, label=false)
plot!(s, tip_loss_function_Pazy("VLM-def", deg2rad(θq), Uq, deg2rad(Λq)).(s), c=:red, label="Interpolated: Λ=$(Λq)°, θ=$(θq)°, U=$Uq m/s")
display(plt)
