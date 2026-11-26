using AeroBeams

# Directory of the static aeroelastic VLM solutions of the Pazy wing: one folder per root angle of attack (aoa<θ>deg, with "p" as decimal separator), each with the spanwise distributions (spanwise_cn_alpha.csv) for all the sweep angles and a range of airspeeds
vlmDir = joinpath(pkgdir(AeroBeams), "test", "referenceData", "PazyVLM")

# Speed of sound used for the compressibility factor [m/s] (that of the airfoil's lift slope in AeroBeams, at sea level)
a = standard_atmosphere(0).a

# Compressibility factor as a function of airspeed
β_of(U) = sqrt(1-(U/a)^2)

# Tip loss types, each one built from a lift slope of the VLM data. The tangent slopes are the change of normal force per change of local incidence of the sections when the root angle is perturbed, with the deformed shape of the wing fixed ("fixed_shape", the slope seen by the perturbations of the aeroelastic problem, as in flutter) or re-deformed at each perturbation ("equilibrium"), from cn_alpha_local_perturbation_per_rad of spanwise_cn_alpha_perturbations.csv. The secant slope is the change of normal force per local incidence at the static equilibrium (cn_alpha_local_per_rad of spanwise_cn_alpha.csv, cn/α_local), which for the swept wings at high airspeed is dominated by the load redistribution caused by the twist
types = (:VLM_tan_fixed, :VLM_sec, :VLM_tan_equil)
perturbationModes = Dict(:VLM_tan_fixed => "fixed_shape", :VLM_tan_equil => "equilibrium")

# Lift slope [1/rad] along the span, set to NaN where it is not valid: the local incidence does not change with the perturbation (ill-conditioned), or the slope is out of any plausible range. These entries are set further below
function plausible_slope(slope,valid)
    slope = Float64.(slope)
    slope[.!(valid .== true) .| .!isfinite.(slope) .| (slope .< 0.05*2π) .| (slope .> 1.3*2π)] .= NaN
    return slope
end

# Read the lift slope along the span and the tip out-of-plane displacement of the VLM data of all conditions (sweep angle [deg], root angle of attack [deg], airspeed [m/s]), for each tip loss type. The conditions in which the VLM does not converge to an equilibrium have no valid lift slope: they are set further below
slopes = Dict(type => Dict{NTuple{3,Float64},Vector{Float64}}() for type in types)
tipOOP = Dict{NTuple{3,Float64},Float64}()
invalid = NTuple{3,Float64}[]
eta = Float64[]
for aoaDir in filter(d -> occursin(r"^aoa[\dp]+deg$", d), readdir(vlmDir))
    # Read data and set the function to get its columns
    data,header = readdlm(joinpath(vlmDir, aoaDir, "spanwise_cn_alpha.csv"), ',', header=true)
    column(name) = data[:, findfirst(==(name), vec(header))]
    pertData,pertHeader = readdlm(joinpath(vlmDir, aoaDir, "spanwise_cn_alpha_perturbations.csv"), ',', header=true)
    pertColumn(name) = pertData[:, findfirst(==(name), vec(pertHeader))]
    θ = only(unique(round.(column("alpha_deg"), digits=6)))
    # Data of each sweep angle and airspeed (the cases are numbered across both)
    for case in unique(column("case"))
        rows = column("case") .== case
        Λ = only(unique(round.(column("sweep_deg")[rows], digits=6)))
        U = only(unique(round.(column("U_m_s")[rows], digits=6)))
        if !all(column("central_converged")[rows] .== true) || !all(column("cn_alpha_local_valid")[rows] .== true)
            println("The VLM did not converge, or has no valid lift slope, at Λ=$(Λ)°, θ=$(θ)°, U=$U m/s ($aoaDir)")
            push!(invalid, (Λ,θ,U))
            continue
        end
        # The span stations are the same for all conditions
        etaCase = column("eta")[rows]
        isempty(eta) ? append!(eta, etaCase) : @assert etaCase ≈ eta "the span stations of $aoaDir (case $case) differ from the others"
        tipOOP[(Λ,θ,U)] = only(unique(column("tip_OOP_m")[rows]))
        slopes[:VLM_sec][(Λ,θ,U)] = plausible_slope(column("cn_alpha_local_per_rad")[rows], column("cn_alpha_local_valid")[rows])
        for type in (:VLM_tan_fixed, :VLM_tan_equil)
            pertRows = (pertColumn("case") .== case) .& (pertColumn("mode") .== perturbationModes[type])
            @assert pertColumn("eta")[pertRows] ≈ etaCase "the span stations of the perturbations of $aoaDir (case $case) differ from those of the equilibrium"
            slopes[type][(Λ,θ,U)] = plausible_slope(pertColumn("cn_alpha_local_perturbation_per_rad")[pertRows], pertColumn("valid")[pertRows])
        end
    end
end
ΛRange = sort(unique(vcat(getindex.(keys(tipOOP), 1), getindex.(invalid, 1))))
θData = sort(unique(vcat(getindex.(keys(tipOOP), 2), getindex.(invalid, 2))))
URange = sort(unique(vcat(getindex.(keys(tipOOP), 3), getindex.(invalid, 3))))

# The VLM data has a condition with root angle of attack and airspeed close to zero (the smallest ones), in which the wing is not deformed: it holds the lift slope of the undeformed wing, and is the anchor of the table. At zero root pitch the wing of AeroBeams is undeformed (without gravity, as in the flutter sweeps at θ=0), and keeps this tip loss function at any airspeed, so the angle of the anchor is set to zero in the table. The VLM data of the same angle at higher airspeeds is not used: it is the linear aeroelastic response of the VLM, which does not apply to the undeformed wing of AeroBeams and is ill-conditioned (even with negative slopes) when the VLM wing is close to divergence
θAnchor, UAnchor = first(θData), first(URange)
@assert θAnchor <= 1e-2 && UAnchor <= 0.1 "the VLM data must have a condition with root angle of attack and airspeed close to zero, as the anchor of the table"
for Λ in ΛRange
    println("Anchor of Λ=$(Λ)°: θ=$(θAnchor)°, U=$(UAnchor) m/s, tip OOP displacement of $(tipOOP[(Λ,θAnchor,UAnchor)]) m")
end
θRange = vcat(0.0, θData[2:end])

# Ranges of the conditions, which must form a complete grid
@assert eta ≈ range(0, 1, length=length(eta)) "the span stations must be evenly spaced from the root to the tip"
println("VLM data at Λ = $ΛRange deg, θ = $θRange deg, U = $URange m/s, and $(length(eta)) span stations")

# Tip loss function of each condition of a tip loss type: its lift slope, completed over the grid of conditions as described below, normalized by that of the airfoil, 2π/β. No sweep factor is needed: the VLM coefficient is based on the chord of the AeroBeams model and on the spanwise width projected on the unswept direction (Δℓ cosΛ), which matches the sweep correction (cosΛ) of the loads of the strips of AeroBeams. Compare the static equilibria of AeroBeams and of the VLM (tip_OOP_m and tip_twist_deg in the data) to check it
function tip_loss_conditions(type,slopes)

    # Conditions in which the VLM does not converge (the wing is far beyond the range of airspeeds and angles of the experiments): use the lift slope of the closest airspeed with data at the same sweep and root angle, as for the airspeeds outside of the range of the data
    for (Λ,θ,U) in invalid
        UValid = [k[3] for k in keys(tipOOP) if k[1] == Λ && k[2] == θ]
        UClosest = UValid[argmin(abs.(UValid .- U))]
        slopes[(Λ,θ,U)] = copy(slopes[(Λ,θ,UClosest)])
        println("Using the lift slope of $type at U=$UClosest m/s for Λ=$(Λ)°, θ=$(θ)°, U=$U m/s")
    end

    # The undeformed wing (zero root angle) has the lift slope of the anchor, at all airspeeds (apart from the compressibility factor)
    @assert all(all(isfinite, slopes[(Λ,θAnchor,UAnchor)]) for Λ in ΛRange) "the anchor of the table must have a valid lift slope"
    for Λ in ΛRange, U in URange
        slopes[(Λ,0.0,U)] = copy(slopes[(Λ,θAnchor,UAnchor)])
    end

    # Entries of the lift slope that are not valid (NaN): those of the undeformed wing, which the slope is close to elsewhere
    nReplaced = 0
    for (key,slope) in slopes
        bad = isnan.(slope)
        slope[bad] .= slopes[(key[1],0.0,key[3])][bad]
        nReplaced += count(bad)
    end
    println("Entries of the lift slope of $type replaced by those of the undeformed wing: $nReplaced")

    # Ranges of the conditions, which must form a complete grid
    @assert all(haskey(slopes, (Λ,θ,U)) for Λ in ΛRange, θ in θRange, U in URange) "the VLM data does not cover all combinations of sweep angle, root angle of attack and airspeed"

    return Dict((Λ,θ,U) => slopes[(Λ,θ,U)] / (2π/β_of(U)) for Λ in ΛRange, θ in θRange, U in URange)

end

# Save the tip loss function of each type at the span stations, one row per condition
conditionsOf = Dict(type => tip_loss_conditions(type,slopes[type]) for type in types)
for type in types
    table = [[Λ θ U transpose(conditionsOf[type][(Λ,θ,U)])] for Λ in ΛRange for θ in θRange for U in URange]
    header = ["sweep_deg" "aoa_deg" "U_m_s" ["s=$(k-1)/$(length(eta)-1)" for k in eachindex(eta)]...]
    open(joinpath(pkgdir(AeroBeams), "test", "referenceData", "sweptPazy", "tipLossFunction_$(type).csv"), "w") do io
        writedlm(io, header, ',')
        writedlm(io, vcat(table...), ',')
    end
end

# Plot, for each type, a grid condition (the data points are exactly reproduced) and a query between grid points, as used by the models. Restart Julia if the data was already loaded
Λq, θq, Uq = 15, 6, 65
Λg, θg, Ug = 20, 5, 60
s = 0:0.01:1
plts = map(types) do type
    plt = plot(xlims=[0,1], ylims=[0,1.5], title=string(type), xlabel="Normalized span", ylabel="Normalized lift slope, \$c_{n_\\alpha}/(2\\pi/\\beta)\$")
    scatter!(eta, conditionsOf[type][(Float64(Λg),Float64(θg),Float64(Ug))], c=:black, label="VLM data: Λ=$(Λg)°, θ=$(θg)°, U=$Ug m/s")
    plot!(s, tip_loss_function_Pazy(type, deg2rad(θg), Ug, deg2rad(Λg)).(s), c=:black, label=false)
    plot!(s, tip_loss_function_Pazy(type, deg2rad(θq), Uq, deg2rad(Λq)).(s), c=:red, label="Interpolated: Λ=$(Λq)°, θ=$(θq)°, U=$Uq m/s")
end
display(plot(plts..., layout=(1,3), size=(1500,450)))
