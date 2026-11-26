# Verify the strip-load correction actually used by each solved native element.
function verify_pazy_vlm_def(model,airspeed,pitchAngle,sweepAngle)
    wing = model.beams[1]
    nodes = wing.discretization.normalizedNodalPositions
    profile = AeroBeams.tip_loss_function_Pazy("VLM-def",pitchAngle,airspeed,sweepAngle)
    for (i,element) in enumerate(wing.elements)
        distribution = element.aero.loadsDistribution
        midpoint = (nodes[i] + nodes[i+1]) / 2
        @assert isapprox(distribution.Uref,airspeed) "Wrong airspeed for the native Pazy tip-loss profile"
        @assert isapprox(distribution.ϖMid,profile(midpoint)) "Wrong VLM-def factor at Pazy element $i"
    end
end
