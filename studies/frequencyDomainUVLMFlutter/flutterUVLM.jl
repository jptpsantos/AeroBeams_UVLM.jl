using AeroBeams, Printf

# Flutter of the undeformed Pazy wing (zero root pitch angle) with the aerodynamic forces of the frequency-domain unsteady vortex-lattice method (UVLM) in the basis of the vibration modes of AeroBeams: a check of the strip theory of AeroBeams (with the tip loss function), as the three-dimensional unsteady potential flow solution includes the effects of the finite wing on the unsteady loads (the lag of the circulatory loads in particular). The vortex-lattice method was validated against the steady lift slope of the VLM data of the Pazy wing (test/referenceData/PazyVLM) and against Theodorsen's solution for a wing of very high aspect ratio
include(joinpath(@__DIR__, "modalData.jl"))
include(joinpath(@__DIR__, "unsteadyVLM.jl"))

# Sweep angle [deg] (0, 10 or 20), tip mass configuration (:baseline, :LE or :TE, for the swept wings), and flag for the wingtip rod to be modeled as a beam (otherwise its inertia is lumped, with the ad hoc correction of its torsional inertia)
sweepAngleDeg = 0
tipMassConfig = :baseline
beamWingtipRod = false

# Vortex lattice: chordwise and spanwise panels, and wake length in chords (the flutter speed is converged for 10 chords or more)
M = 8
N = 30
wakeChords = 20.0

# Airspeed range [m/s] and air density [kg/m³]
URange = collect(40.0:0.5:110.0)
ρ = 1.225

# Reduced frequencies (k = ωb/U) at which the aerodynamic forces are computed (they are interpolated linearly)
kRange = vcat([0.002, 0.004, 0.007], collect(0.01:0.01:0.1), collect(0.12:0.02:0.6), [0.7, 0.8, 1.0, 1.3, 1.6, 2.0])

# Flag to also compute the flutter with the strip theory of AeroBeams, with and without the finite-wing lag correction
compareWithAeroBeams = true
airfoil = deepcopy(FlatPlate)

## Modes and vortex lattice
Λ = deg2rad(sweepAngleDeg)
data = pazy_modal_data(sweepAngle=Λ, tipMassConfig=tipMassConfig, beamWingtipRod=beamWingtipRod)
c, xs, L = data.chord, data.chord*data.referenceAxis, data.length
b = c/2

# Modes with aerodynamic effect (the in-plane modes do not move the wing in the lift direction nor twist it)
modes = filter(m -> maximum(abs, m.h) > 5e-6 || maximum(abs, m.θ) > 1e-4, data.modes)
println("Frequencies of the modes [Hz]: ", round.([m.frequency for m in modes], digits=1))
lattice = create_vortex_lattice(c=c, Lp=L*cos(Λ), Λ=Λ, xs0=xs, M=M, N=N, wakeChords=wakeChords)
println("Steady lift slope of the planform, CL/(2πα): ", round(steady_lift_slope(lattice), digits=4))

## Aerodynamic forces and flutter
fields = mode_fields(modes, data.s)
Qtable = [generalized_forces(lattice, fields, L, k) for k in kRange]
function Qinterp(k)
    k <= kRange[1] && return Qtable[1]
    k >= kRange[end] && return Qtable[end]
    i = searchsortedlast(kRange, k)
    t = (k - kRange[i])/(kRange[i+1] - kRange[i])
    return (1-t)*Qtable[i] + t*Qtable[i+1]
end
roots = pk_flutter([m.M for m in modes], [m.K for m in modes], Qinterp, b, ρ, URange)
println("\nUVLM (M=$M, N=$N, wake of $wakeChords chords): branches of the p-k solution that cross the stability limit")
for r in 1:min(6, length(modes))
    onset = stability_crossings(URange, roots[:,r])
    offset = stability_crossings(URange, roots[:,r], upward=false)
    isempty(onset) && isempty(offset) && continue
    iU = argmin(abs.(URange .- (isempty(onset) ? offset[1] : onset[1])))
    @printf("  branch of mode %d (%.1f Hz at zero airspeed): onset %s m/s, offset %s m/s, frequency at the crossing %.1f Hz\n", r, imag(roots[1,r])/2π, string(round.(onset, digits=1)), string(round.(offset, digits=1)), imag(roots[iU,r])/2π)
end
println("\nFrequency [Hz] and growth rate [1/s] of the first four branches")
for U in 50.0:10.0:100.0
    iU = argmin(abs.(URange .- U))
    @printf("  U = %5.1f m/s:", U)
    for r in 1:4
        @printf("  f=%6.2f σ=%7.2f |", imag(roots[iU,r])/2π, real(roots[iU,r]))
    end
    println()
end

## Strip theory of AeroBeams (eigenproblem at zero pitch angle)
if compareWithAeroBeams
    nElem,_,_,_,_ = geometrical_properties_Pazy(Λ=Λ)
    for finiteWingLag in (false, true)
        eigensols = Vector{EigenSolution}(undef, length(URange))
        x0 = zeros(0)
        for (j,U) in enumerate(URange)
            kwargs = Λ == 0 ? (;) : (; tipMassConfig=tipMassConfig, beamWingtipRod=beamWingtipRod, tipRodInertiaCorrection=!beamWingtipRod)
            model,_ = create_Pazy(aeroSolver=Indicial(finiteWingLag=finiteWingLag), airfoil=deepcopy(airfoil), sweepAngle=Λ, pitchAngle=0, airspeed=U; kwargs...)
            problem = create_EigenProblem(model=model, eigenmethod=:krylovkit, nev=Λ == 0 ? 45 : 40, x0=x0)
            solve!(problem)
            eigensols[j] = problem.eigensolution
            x0 = problem.systemArrays.x
        end
        freqs,damps,_ = mode_tracking(eigensols)
        rads2Hz!(freqs)
        freqsOsc,dampsOsc,_ = filter_tracked_modes(freqs, damps; fmin=1e-1, maxTracks=5)
        onset,onsetFreq,_,_,offset,_,_,_ = find_flutter(URange, freqsOsc, dampsOsc, zeros(length(URange)))
        println("AeroBeams (flat plate, VLM, finiteWingLag=$finiteWingLag): onset ", round.(onset, digits=1), " m/s at ", round.(onsetFreq, digits=1), " Hz, offset ", round.(offset, digits=1), " m/s")
    end
end
