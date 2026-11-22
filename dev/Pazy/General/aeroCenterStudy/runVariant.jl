# Runs the calculation part (up to "## Plots") of a flutter or static script of the Pazy wing with beam rod, with the variant of the aerodynamic model selected by environment variables (see variants.jl), and prints the errors relative to the experiments:
#   julia -t 1 --project=<AeroBeams> runVariant.jl <script> <tag>
# <script>: s0 (S0/BeamRod/S0BeamRodHumpFlutterRegion.jl), s10 (S10/BeamRod/S10BeamRodFlutterOnsetPitchRange.jl), s20 (S20/BeamRod/S20BeamRodFlutterOnsetPitchRange.jl) or s0static (S0/BeamRod/S0BeamRodPitchRange.jl). The raw metrics of the run are saved in outputs/<script>_<tag>.txt
# Example (airfoil of the data of Timmer with the tangent lift slope in the eigenproblem): AF=timmerfull TLMODE=swap julia -t 1 --project=. dev/Pazy/General/aeroCenterStudy/runVariant.jl s10 timmer_swap
using Printf, DelimitedFiles
using AeroBeams

const SCRIPTS = Dict(
    "s0"       => "dev/Pazy/S0/BeamRod/S0BeamRodHumpFlutterRegion.jl",
    "s10"      => "dev/Pazy/S10/BeamRod/S10BeamRodFlutterOnsetPitchRange.jl",
    "s20"      => "dev/Pazy/S20/BeamRod/S20BeamRodFlutterOnsetPitchRange.jl",
    "s0static" => "dev/Pazy/S0/BeamRod/S0BeamRodPitchRange.jl",
)
key, tag = ARGS[1], ARGS[2]
OVERRIDES = Pair{String,String}[]
include(joinpath(@__DIR__, "variants.jl"))

# Calculation part of the script, with the replacements of the variant
# (the scripts are those of the repository containing this folder, which may differ from the loaded package, e.g., a copy of it with the prototype of aeroCenterShift.patch)
script = joinpath(abspath(joinpath(@__DIR__, "..", "..", "..", "..")), SCRIPTS[key])
txt = read(script, String)
txt = txt[1:first(findfirst("## Plots", txt))-1]
txt = foldl(OVERRIDES; init=txt) do t, (a,b)
    occursin(a, t) || error("replacement pattern not found in $(SCRIPTS[key]): $a")
    replace(t, a => b)
end
t0 = time()
include_string(Main, txt, script)
elapsed = time() - t0

## Errors relative to the experiments
interp(x, y, xq) = (xq < first(x) || xq > last(x)) ? NaN : begin
    j = min(searchsortedlast(x, xq), length(x)-1); y[j] + (y[j+1]-y[j])*(xq-x[j])/(x[j+1]-x[j])
end
rms(v) = sqrt(sum(abs2, v)/length(v))
mkpath(joinpath(@__DIR__, "outputs"))
io = open(joinpath(@__DIR__, "outputs", "$(key)_$(tag).txt"), "w")
report(args...) = (println(args...); println(io, args...))
report("# $(SCRIPTS[key]) [$tag]: AF=$(AF), FWL=$(get(ENV,"FWL","default")), TLMODE=$tlmode, XS=$useShift, WS=$WSCALE, WN=$useNastran ($(round(elapsed,digits=1)) s)")
if key == "s0"
    # Hump flutter region: first onset and offset of the hump mode (mode 3) below 100 m/s, vs. the onset of the tests when the airspeed is increased (lower edge) and decreased (upper edge)
    report("θ [deg] | onset  offset | exp. onset up / down | errors (m/s)")
    errUp, errDown = Float64[], Float64[]
    for (i,θ) in enumerate(θRange)
        on = findfirst(m -> m == 3, flutterOnsetModes[i]); off = findfirst(m -> m == 3, flutterOffsetModes[i])
        Uon = isnothing(on) || flutterOnsetSpeeds[i][on] >= 100 ? NaN : flutterOnsetSpeeds[i][on]
        Uoff = isnothing(off) ? NaN : flutterOffsetSpeeds[i][off]
        iUp, iDown = findfirst(≈(θ), rootPitchVelUp), findfirst(≈(θ), rootPitchVelDown)
        eUp = isnothing(iUp) ? NaN : Uon - flutterOnsetVelUp[iUp]
        eDown = isnothing(iDown) ? NaN : Uoff - flutterOnsetVelDown[iDown]
        push!(errUp, eUp); push!(errDown, eDown)
        fmt(e) = isnan(e) ? "    -" : @sprintf("%+5.1f", e)
        report(@sprintf("%6.1f  | %6.1f %6.1f | %7s %7s | %s %s", θ, Uon, Uoff, isnothing(iUp) ? "-" : string(flutterOnsetVelUp[iUp]), isnothing(iDown) ? "-" : string(flutterOnsetVelDown[iDown]), fmt(eUp), fmt(eDown)))
    end
    report(@sprintf("mean error lower edge %+.2f, upper edge %+.2f m/s", sum(filter(!isnan, errUp))/count(!isnan, errUp), sum(filter(!isnan, errDown))/count(!isnan, errDown)))
elseif key == "s10" || key == "s20"
    # First flutter onset vs. root pitch angle, per tip mass configuration
    allErr = Float64[]
    for (c,config) in enumerate(tipMassConfigs)
        ref = config == :LE ? Uf_vs_aoa_LE : Uf_vs_aoa_TE
        refFreq = config == :LE ? omegaf_vs_aoa_LE : omegaf_vs_aoa_TE
        errs, errsf = Float64[], Float64[]
        line, linef = "", ""
        for (θ,Uexp,fexp) in zip(ref[1,:], ref[2,:], refFreq[2,:])
            i = findfirst(x -> isapprox(rad2deg(x), θ; atol=1e-6), θRange)
            isnothing(i) && continue
            push!(errs, flutterOnsetSpeed[c,i] - Uexp); push!(errsf, rads2Hz(flutterOnsetFreq[c,i]) - fexp)
            line *= @sprintf(" %g:%+.1f", θ, errs[end]); linef *= @sprintf(" %g:%+.1f", θ, errsf[end])
        end
        append!(allErr, errs)
        report("$(uppercase(key)) $config: flutter speed error [m/s] at θ [deg]:$line | mean $(@sprintf("%+.2f", sum(errs)/length(errs))) | frequency error [Hz]:$linef")
    end
    report(@sprintf("mean absolute error of the flutter speeds: %.2f m/s", sum(abs, allErr)/length(allErr)))
elseif key == "s0static"
    # Tip OOP displacement vs. dynamic pressure (test 1) and tip twist vs. airspeed of the experiments
    ρ = problem[1].model.atmosphere.ρ
    for (i,θ) in enumerate(θRange)
        θd = round(Int, rad2deg(θ))
        ref = getfield(Main, Symbol("tip_u3Vsq_rootPitch$(θd)_ExpTest1"))
        Uref = sqrt.(2 .* ref[1,:] ./ ρ)
        zab = [interp(collect(URange), tip_OOP[i,:], u) for u in Uref]
        ok = .!isnan.(zab) .& (Uref .>= 20)
        line = @sprintf("θ=%d°: tip OOP rms %.2f %%L (mean relative error %+.1f %%)", θd, rms((zab[ok] .- ref[2,:][ok]) ./ L .* 100), 100*sum((zab[ok] .- ref[2,:][ok]) ./ ref[2,:][ok])/count(ok))
        if θd in (5,7)
            tw = getfield(Main, Symbol("tip_thetaVsU_rootPitch$(θd)_Exp"))
            twab = [interp(collect(URange), tip_twist[i,:], u) for u in tw[1,:]]
            ok2 = .!isnan.(twab) .& (tw[1,:] .>= 20)
            line *= @sprintf(" | twist rms %.2f deg (mean relative error %+.1f %%)", rms(twab[ok2] .- tw[2,:][ok2]), 100*sum((twab[ok2] .- tw[2,:][ok2]) ./ tw[2,:][ok2])/count(ok2))
        end
        report(line)
    end
end
close(io)
