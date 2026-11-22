using AeroBeams

# Root pitch angle [rad] and airspeed [m/s] ranges
θRange = deg2rad.([3,5,7])
URange = vcat(0:1:60)

# Aerodynamic solver
aeroSolver = Indicial()

# Mesh (defines number of elements)
mesh = :nominal

# Flag to use AePW3 matrices for sectional stiffness
useAePW3Matrices = false

# Fixed geometrical and discretization properties
nElem,L,chord,referenceAxis,nodalPositions = geometrical_properties_Pazy(; mesh=mesh, useAePW3Matrices=useAePW3Matrices)

# Tip mass (for test 1, 0.01 kg, 40 mm behind the trailing-edge)
tipMass = 0.01
ηtipMass = [0; -chord*(1-referenceAxis)-0.04; 0]

# Initialize outputs
problem = Array{SteadyProblem}(undef,length(θRange),length(URange))
tip_OOP = Array{Float64}(undef,length(θRange),length(URange))
tip_twist = Array{Float64}(undef,length(θRange),length(URange))

# Sweep root angle
for (i,θ) in enumerate(θRange)
    # Initialize guess solution
    x0 = zeros(0)
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        # Display progress
        println("Solving for θ = $(round(rad2deg(θ))) deg, U = $U m/s")
        # Update model
        PazyWingPitchRange,_ = create_Pazy(aeroSolver=aeroSolver,pitchAngle=θ,airspeed=U,tipMasses=[tipMass],ηtipMasses=[ηtipMass],tipMassesInertias=[zeros(3,3)],mesh=mesh,useAePW3Matrices=useAePW3Matrices)
        # Create and solve problem
        problem[i,j] = create_SteadyProblem(model=PazyWingPitchRange,x0=x0)
        solve!(problem[i,j])
        # Get tip twist, AoA, IP and OOP displacement at beam reference line
        states = problem[i,j].loadstepHistory.states
        tipElemFinalAero = problem[i,j].loadstepHistory.aeroVariables[end][nElem]
        tipNodeFinalState = states[end][nElem].nodal
        tip_twist[i,j] = twist_angle(tipNodeFinalState.p2_b)
        tip_OOP[i,j] = -tipNodeFinalState.u2[1]
        # Set initial guess solution for next iteration
        x0 = problem[i,j].systemArrays.x
    end
end

# Load reference data
tip_u3Vsq_rootPitch3_ExpTest1 = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Pazy/tip_u3Vsq_rootPitch3_ExpTest1.txt")
tip_u3Vsq_rootPitch5_ExpTest1 = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Pazy/tip_u3Vsq_rootPitch5_ExpTest1.txt")
tip_u3Vsq_rootPitch7_ExpTest1 = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Pazy/tip_u3Vsq_rootPitch7_ExpTest1.txt")
tip_thetaVsU_rootPitch5_Exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Pazy/tip_thetaVsU_rootPitch5_Exp.txt")
tip_thetaVsU_rootPitch7_Exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Pazy/tip_thetaVsU_rootPitch7_Exp.txt")

## Plots
# Set paths
relPath = "/dev/Pazy/S0/InertiaRod/outputs/PazyWingPitchRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow,length(θRange), categorical=true)
ts = 10
fs = 16
lfs = 9
lw = 2
ms = 3
msw = 0
gr()

# Dynamic pressure range
ρ = problem[1].model.atmosphere.ρ
qRange = 1/2*ρ*URange.^2

# Tip OOP displacement vs. dynamic pressure for root several pitch angles 
plt1 = plot(xlabel="Dynamic pressure [Pa]", ylabel="Tip OOP disp. [% semispan]", xlims=[0,2200], ylims=[0,55], tickfont=font(ts), guidefont=font(fs), legend=:topleft, legendfontsize=lfs)
plot!([NaN], [NaN], c=:black, lw=lw, ls=:solid, label="AeroBeams")
scatter!([NaN], [NaN], c=:black, shape=:diamond, ms=ms, label="Exp. - Avin et al. (2022)")
for (i,θ) in enumerate(θRange)
    plot!(qRange, tip_OOP[i,:]/L*100, c=colors[i], lw=lw, ls=:solid, label=string("\$\\alpha_r=",round(Int,rad2deg(θ)),"^\\circ\$"))
    if i==1
        scatter!(tip_u3Vsq_rootPitch3_ExpTest1[1,:], tip_u3Vsq_rootPitch3_ExpTest1[2,:]*100/L, mc=colors[i], shape=:diamond, ms=ms, msw=msw, label=false)
    elseif i==2
        scatter!(tip_u3Vsq_rootPitch5_ExpTest1[1,:], tip_u3Vsq_rootPitch5_ExpTest1[2,:]*100/L, mc=colors[i], shape=:diamond, ms=ms, msw=msw, label=false)
    elseif i==3
        scatter!(tip_u3Vsq_rootPitch7_ExpTest1[1,:], tip_u3Vsq_rootPitch7_ExpTest1[2,:]*100/L, mc=colors[i], shape=:diamond, ms=ms, msw=msw, label=false)
    end
end
display(plt1)
savefig(string(absPath,"/PazyWingPitchRange_tipOOP.pdf"))

# Tip twist vs. airspeed for root several pitch angles 
plt2 = plot(xlabel="Airspeed [m/s]", ylabel="Tip twist [deg]", xlims=[0,60], tickfont=font(ts), guidefont=font(fs))
for (i,θ) in enumerate(θRange)
    plot!(URange, tip_twist[i,:], c=colors[i], lw=lw, label=false)
    if i==2
        scatter!(tip_thetaVsU_rootPitch5_Exp[1,:], tip_thetaVsU_rootPitch5_Exp[2,:], mc=colors[i], shape=:diamond, ms=ms, msw=msw, label=false)
    elseif i==3
        scatter!(tip_thetaVsU_rootPitch7_Exp[1,:], tip_thetaVsU_rootPitch7_Exp[2,:], mc=colors[i], shape=:diamond, ms=ms, msw=msw, label=false)
    end
end
display(plt2)
savefig(string(absPath,"/PazyWingPitchRange_tipTwist.pdf"))

println("Finished PazyWingPitchRange.jl")