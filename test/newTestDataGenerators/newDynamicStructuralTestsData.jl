# Dynamic analysis of the axial vibration of a beam under a traction force applied suddenly
include("../examples/dynamicStructural/axialTractionCantilever.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicStructural/axialTractionCantilever"))
writedlm("test/newTestDataGenerators/dynamicStructural/axialTractionCantilever/u1_08.txt", test_data_matrix(u1_08))
writedlm("test/newTestDataGenerators/dynamicStructural/axialTractionCantilever/u1_10.txt", test_data_matrix(u1_10))

# Dynamic analysis of the free response of a beam clamped at both ends and subjected to an initial displacement profile
include("../examples/dynamicStructural/biclampedBeam.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicStructural/biclampedBeam"))
writedlm("test/newTestDataGenerators/dynamicStructural/biclampedBeam/u3_mid.txt", test_data_matrix(u3_mid))
writedlm("test/newTestDataGenerators/dynamicStructural/biclampedBeam/V3_mid.txt", test_data_matrix(V3_mid))
writedlm("test/newTestDataGenerators/dynamicStructural/biclampedBeam/Vdot3_mid.txt", test_data_matrix(Vdot3_mid))
writedlm("test/newTestDataGenerators/dynamicStructural/biclampedBeam/theta2_quarter.txt", test_data_matrix(θ2_quarter))
writedlm("test/newTestDataGenerators/dynamicStructural/biclampedBeam/Omega2_quarter.txt", test_data_matrix(Ω2_quarter))
writedlm("test/newTestDataGenerators/dynamicStructural/biclampedBeam/Omegadot2_quarter.txt", test_data_matrix(Ωdot2_quarter))

# Dynamic analysis of the free decay of a cantilever with stiffness-proportional structural damping
include("../examples/dynamicStructural/cantileverStructuralDampingFreeDecay.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicStructural/cantileverStructuralDampingFreeDecay"))
writedlm("test/newTestDataGenerators/dynamicStructural/cantileverStructuralDampingFreeDecay/u3_tip.txt", test_data_matrix(u3Tip))

# Dynamic analysis of a slanted double pendulum released from rest
include("../examples/dynamicStructural/ChengSlantedDoublePendulum.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicStructural/ChengSlantedDoublePendulum"))
writedlm("test/newTestDataGenerators/dynamicStructural/ChengSlantedDoublePendulum/u1_tip.txt", test_data_matrix(u1_tip))
writedlm("test/newTestDataGenerators/dynamicStructural/ChengSlantedDoublePendulum/u2_tip.txt", test_data_matrix(u2_tip))
writedlm("test/newTestDataGenerators/dynamicStructural/ChengSlantedDoublePendulum/u3_tip.txt", test_data_matrix(u3_tip))

# Dynamic analysis of a composite cantilever beam under a tip sinusoidal load
include("../examples/dynamicStructural/compositeCantilever.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicStructural/compositeCantilever"))
writedlm("test/newTestDataGenerators/dynamicStructural/compositeCantilever/u1_tip.txt", test_data_matrix(u1_tip))
writedlm("test/newTestDataGenerators/dynamicStructural/compositeCantilever/u2_tip.txt", test_data_matrix(u2_tip))
writedlm("test/newTestDataGenerators/dynamicStructural/compositeCantilever/u3_tip.txt", test_data_matrix(u3_tip))
writedlm("test/newTestDataGenerators/dynamicStructural/compositeCantilever/p1_tip.txt", test_data_matrix(p1_tip))
writedlm("test/newTestDataGenerators/dynamicStructural/compositeCantilever/p2_tip.txt", test_data_matrix(p2_tip))
writedlm("test/newTestDataGenerators/dynamicStructural/compositeCantilever/p3_tip.txt", test_data_matrix(p3_tip))
writedlm("test/newTestDataGenerators/dynamicStructural/compositeCantilever/F1_root.txt", test_data_matrix(F1_root))
writedlm("test/newTestDataGenerators/dynamicStructural/compositeCantilever/F2_root.txt", test_data_matrix(F2_root))
writedlm("test/newTestDataGenerators/dynamicStructural/compositeCantilever/F3_root.txt", test_data_matrix(F3_root))
writedlm("test/newTestDataGenerators/dynamicStructural/compositeCantilever/M1_root.txt", test_data_matrix(M1_root))
writedlm("test/newTestDataGenerators/dynamicStructural/compositeCantilever/M2_root.txt", test_data_matrix(M2_root))
writedlm("test/newTestDataGenerators/dynamicStructural/compositeCantilever/M3_root.txt", test_data_matrix(M3_root))

# Dynamic analysis of a curved cantilever subjected to a tip follower force
include("../examples/dynamicStructural/curvedCantileverDynamicFollower.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicStructural/curvedCantileverDynamicFollower"))
writedlm("test/newTestDataGenerators/dynamicStructural/curvedCantileverDynamicFollower/u1_tip.txt", test_data_matrix(u1_tip))
writedlm("test/newTestDataGenerators/dynamicStructural/curvedCantileverDynamicFollower/u2_tip.txt", test_data_matrix(u2_tip))
writedlm("test/newTestDataGenerators/dynamicStructural/curvedCantileverDynamicFollower/u3_tip.txt", test_data_matrix(u3_tip))
writedlm("test/newTestDataGenerators/dynamicStructural/curvedCantileverDynamicFollower/F1_root.txt", test_data_matrix(F1_root))
writedlm("test/newTestDataGenerators/dynamicStructural/curvedCantileverDynamicFollower/F2_root.txt", test_data_matrix(F2_root))
writedlm("test/newTestDataGenerators/dynamicStructural/curvedCantileverDynamicFollower/F3_root.txt", test_data_matrix(F3_root))
writedlm("test/newTestDataGenerators/dynamicStructural/curvedCantileverDynamicFollower/M1_root.txt", test_data_matrix(M1_root))
writedlm("test/newTestDataGenerators/dynamicStructural/curvedCantileverDynamicFollower/M2_root.txt", test_data_matrix(M2_root))
writedlm("test/newTestDataGenerators/dynamicStructural/curvedCantileverDynamicFollower/M3_root.txt", test_data_matrix(M3_root))

# Dynamic analysis of a double pendulum released from rest
include("../examples/dynamicStructural/doublePendulum.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicStructural/doublePendulum"))
writedlm("test/newTestDataGenerators/dynamicStructural/doublePendulum/u1_hinge.txt", test_data_matrix(u1_hinge))
writedlm("test/newTestDataGenerators/dynamicStructural/doublePendulum/u3_hinge.txt", test_data_matrix(u3_hinge))
writedlm("test/newTestDataGenerators/dynamicStructural/doublePendulum/u1_tip.txt", test_data_matrix(u1_tip))
writedlm("test/newTestDataGenerators/dynamicStructural/doublePendulum/u3_tip.txt", test_data_matrix(u3_tip))

# Dynamic analysis of a right-angled frame subjected to an out-of-plane force
include("../examples/dynamicStructural/elbowFrame.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicStructural/elbowFrame"))
writedlm("test/newTestDataGenerators/dynamicStructural/elbowFrame/u3_elbow.txt", test_data_matrix(u3_elbow))
writedlm("test/newTestDataGenerators/dynamicStructural/elbowFrame/u3_tip.txt", test_data_matrix(u3_tip))

# Dynamic analysis of a very flexible beam subjected to loads yielding two-dimensional motion
include("../examples/dynamicStructural/flyingFlexibleBeam2D.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicStructural/flyingFlexibleBeam2D"))
writedlm("test/newTestDataGenerators/dynamicStructural/flyingFlexibleBeam2D/u1_tip.txt", test_data_matrix(u1_tip))
writedlm("test/newTestDataGenerators/dynamicStructural/flyingFlexibleBeam2D/u3_tip.txt", test_data_matrix(u3_tip))

# Dynamic analysis of a hinged beam in free flight
include("../examples/dynamicStructural/flyingScissors.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicStructural/flyingScissors"))
writedlm("test/newTestDataGenerators/dynamicStructural/flyingScissors/u1_tipA.txt", test_data_matrix(u1_tipA))
writedlm("test/newTestDataGenerators/dynamicStructural/flyingScissors/u3_tipA.txt", test_data_matrix(u3_tipA))
writedlm("test/newTestDataGenerators/dynamicStructural/flyingScissors/u1_tipB.txt", test_data_matrix(u1_tipB))
writedlm("test/newTestDataGenerators/dynamicStructural/flyingScissors/u3_tipB.txt", test_data_matrix(u3_tipB))
writedlm("test/newTestDataGenerators/dynamicStructural/flyingScissors/u1_hinge.txt", test_data_matrix(u1_hinge))
writedlm("test/newTestDataGenerators/dynamicStructural/flyingScissors/u3_hinge.txt", test_data_matrix(u3_hinge))

# Dynamic analysis of a very flexible beam subjected to loads yielding two-dimensional motion
include("../examples/dynamicStructural/flyingSpaghetti2D.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicStructural/flyingSpaghetti2D"))
writedlm("test/newTestDataGenerators/dynamicStructural/flyingSpaghetti2D/u1_tip.txt", test_data_matrix(u1_tip))
writedlm("test/newTestDataGenerators/dynamicStructural/flyingSpaghetti2D/u3_tip.txt", test_data_matrix(u3_tip))

# Dynamic analysis of a very flexible beam subjected to loads yielding tri-dimensional motion
include("../examples/dynamicStructural/flyingSpaghetti3D.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicStructural/flyingSpaghetti3D"))
writedlm("test/newTestDataGenerators/dynamicStructural/flyingSpaghetti3D/u1_tip.txt", test_data_matrix(u1_tip))
writedlm("test/newTestDataGenerators/dynamicStructural/flyingSpaghetti3D/u2_tip.txt", test_data_matrix(u2_tip))
writedlm("test/newTestDataGenerators/dynamicStructural/flyingSpaghetti3D/u3_tip.txt", test_data_matrix(u3_tip))

# Dynamic analysis of the free response of a beam subjected to initial displacement and velocity profiles
include("../examples/dynamicStructural/initialDispAndVelBeam.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicStructural/initialDispAndVelBeam"))
writedlm("test/newTestDataGenerators/dynamicStructural/initialDispAndVelBeam/u3_quarter.txt", test_data_matrix(u3_quarter))
writedlm("test/newTestDataGenerators/dynamicStructural/initialDispAndVelBeam/V3_quarter.txt", test_data_matrix(V3_quarter))
writedlm("test/newTestDataGenerators/dynamicStructural/initialDispAndVelBeam/Vdot3_quarter.txt", test_data_matrix(Vdot3_quarter))
writedlm("test/newTestDataGenerators/dynamicStructural/initialDispAndVelBeam/theta2_root.txt", test_data_matrix(θ2_root))
writedlm("test/newTestDataGenerators/dynamicStructural/initialDispAndVelBeam/Omega2_mid.txt", test_data_matrix(Ω2_mid))
writedlm("test/newTestDataGenerators/dynamicStructural/initialDispAndVelBeam/Omegadot2_mid.txt", test_data_matrix(Ωdot2_mid))

# Dynamic analysis of the free response of a beam subjected to an initial displacement profile
include("../examples/dynamicStructural/initialDisplacementBeam.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicStructural/initialDisplacementBeam"))
writedlm("test/newTestDataGenerators/dynamicStructural/initialDisplacementBeam/u3_quarter.txt", test_data_matrix(u3_quarter))
writedlm("test/newTestDataGenerators/dynamicStructural/initialDisplacementBeam/V3_quarter.txt", test_data_matrix(V3_quarter))
writedlm("test/newTestDataGenerators/dynamicStructural/initialDisplacementBeam/Vdot3_quarter.txt", test_data_matrix(Vdot3_quarter))
writedlm("test/newTestDataGenerators/dynamicStructural/initialDisplacementBeam/theta2_root.txt", test_data_matrix(θ2_root))
writedlm("test/newTestDataGenerators/dynamicStructural/initialDisplacementBeam/Omega2_mid.txt", test_data_matrix(Ω2_mid))
writedlm("test/newTestDataGenerators/dynamicStructural/initialDisplacementBeam/Omegadot2_mid.txt", test_data_matrix(Ωdot2_mid))

# Dynamic analysis of the free response of a beam subjected to an initial velocity profile
include("../examples/dynamicStructural/initialVelocityBeam.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicStructural/initialVelocityBeam"))
writedlm("test/newTestDataGenerators/dynamicStructural/initialVelocityBeam/u3_quarter.txt", test_data_matrix(u3_quarter))
writedlm("test/newTestDataGenerators/dynamicStructural/initialVelocityBeam/V3_quarter.txt", test_data_matrix(V3_quarter))
writedlm("test/newTestDataGenerators/dynamicStructural/initialVelocityBeam/Vdot3_quarter.txt", test_data_matrix(Vdot3_quarter))
writedlm("test/newTestDataGenerators/dynamicStructural/initialVelocityBeam/theta2_root.txt", test_data_matrix(θ2_root))
writedlm("test/newTestDataGenerators/dynamicStructural/initialVelocityBeam/Omega2_mid.txt", test_data_matrix(Ω2_mid))
writedlm("test/newTestDataGenerators/dynamicStructural/initialVelocityBeam/Omegadot2_mid.txt", test_data_matrix(Ωdot2_mid))

# Dynamic analysis of joined beams under load
include("../examples/dynamicStructural/joinedBeams.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicStructural/joinedBeams"))
writedlm("test/newTestDataGenerators/dynamicStructural/joinedBeams/u1.txt", test_data_matrix(u1))
writedlm("test/newTestDataGenerators/dynamicStructural/joinedBeams/u2.txt", test_data_matrix(u2))
writedlm("test/newTestDataGenerators/dynamicStructural/joinedBeams/u3.txt", test_data_matrix(u3))

# Dynamic analysis of a pinned robot arm driven by a couple moment
include("../examples/dynamicStructural/momentDrivenRobotArm.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicStructural/momentDrivenRobotArm"))
writedlm("test/newTestDataGenerators/dynamicStructural/momentDrivenRobotArm/u1_tip.txt", test_data_matrix(u1_tip))
writedlm("test/newTestDataGenerators/dynamicStructural/momentDrivenRobotArm/u3_tip.txt", test_data_matrix(u3_tip))

# Dynamic analysis of a pendulum released from rest
include("../examples/dynamicStructural/pendulum.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicStructural/pendulum"))
writedlm("test/newTestDataGenerators/dynamicStructural/pendulum/u1_tip.txt", test_data_matrix(u1_tip))
writedlm("test/newTestDataGenerators/dynamicStructural/pendulum/u3_tip.txt", test_data_matrix(u3_tip))

# Dynamic analysis of a clamped beam with root sinusoidal oscillation (ω = 9 Hz ≈ 2nd bending mode)
include("../examples/dynamicStructural/rootExcitationBeam.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicStructural/rootExcitationBeam"))
writedlm("test/newTestDataGenerators/dynamicStructural/rootExcitationBeam/u3b_root.txt", test_data_matrix(u3b_root))
writedlm("test/newTestDataGenerators/dynamicStructural/rootExcitationBeam/u3b_tip.txt", test_data_matrix(u3b_tip))
writedlm("test/newTestDataGenerators/dynamicStructural/rootExcitationBeam/V3_root.txt", test_data_matrix(V3_root))
writedlm("test/newTestDataGenerators/dynamicStructural/rootExcitationBeam/V3_tip.txt", test_data_matrix(V3_tip))

# Dynamic analysis of a rotary shaft with specified rotation
include("../examples/dynamicStructural/rotaryShaft.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicStructural/rotaryShaft"))
writedlm("test/newTestDataGenerators/dynamicStructural/rotaryShaft/pNum.txt", test_data_matrix(pNum))
writedlm("test/newTestDataGenerators/dynamicStructural/rotaryShaft/pdotNum.txt", test_data_matrix(pdotNum))
writedlm("test/newTestDataGenerators/dynamicStructural/rotaryShaft/OmegaNum.txt", test_data_matrix(ΩNum))
writedlm("test/newTestDataGenerators/dynamicStructural/rotaryShaft/OmegadotNum.txt", test_data_matrix(ΩdotNum))
writedlm("test/newTestDataGenerators/dynamicStructural/rotaryShaft/MNum.txt", test_data_matrix(MNum))

# Dynamic analysis of an articulated robot arm driven by specified rotation
include("../examples/dynamicStructural/rotationDrivenArticulatedRobotArm.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicStructural/rotationDrivenArticulatedRobotArm"))
writedlm("test/newTestDataGenerators/dynamicStructural/rotationDrivenArticulatedRobotArm/u1_tip.txt", test_data_matrix(u1_tip))
writedlm("test/newTestDataGenerators/dynamicStructural/rotationDrivenArticulatedRobotArm/u3_tip.txt", test_data_matrix(u3_tip))
writedlm("test/newTestDataGenerators/dynamicStructural/rotationDrivenArticulatedRobotArm/u1_hinge.txt", test_data_matrix(u1_hinge))
writedlm("test/newTestDataGenerators/dynamicStructural/rotationDrivenArticulatedRobotArm/u3_hinge.txt", test_data_matrix(u3_hinge))

# Dynamic analysis of a robot arm driven by specified rotation
include("../examples/dynamicStructural/rotationDrivenRobotArm.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicStructural/rotationDrivenRobotArm"))
writedlm("test/newTestDataGenerators/dynamicStructural/rotationDrivenRobotArm/u1_tip.txt", test_data_matrix(u1_tip))
writedlm("test/newTestDataGenerators/dynamicStructural/rotationDrivenRobotArm/u3_tip.txt", test_data_matrix(u3_tip))

# Dynamic analysis of the spin-up maneuver of a robot arm
include("../examples/dynamicStructural/spinupRobotArm.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicStructural/spinupRobotArm"))
writedlm("test/newTestDataGenerators/dynamicStructural/spinupRobotArm/u1_tip.txt", test_data_matrix(u1_tip))
writedlm("test/newTestDataGenerators/dynamicStructural/spinupRobotArm/u2_tip.txt", test_data_matrix(u2_tip))
writedlm("test/newTestDataGenerators/dynamicStructural/spinupRobotArm/theta3_root.txt", test_data_matrix(θ3_root))

# Dynamic analysis of a cantilever with linearly increasing tip moment
include("../examples/dynamicStructural/tipMomentCantileverDynamic.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicStructural/tipMomentCantileverDynamic"))
writedlm("test/newTestDataGenerators/dynamicStructural/tipMomentCantileverDynamic/kappa2.txt", test_data_matrix(κ2))

# Dynamic analysis of a pendulum with tip mass
include("../examples/dynamicStructural/tipPendulum.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicStructural/tipPendulum"))
writedlm("test/newTestDataGenerators/dynamicStructural/tipPendulum/u1_tip.txt", test_data_matrix(u1_tip))
writedlm("test/newTestDataGenerators/dynamicStructural/tipPendulum/u3_tip.txt", test_data_matrix(u3_tip))

# Dynamic analysis of a cantilever with tip sinusoidal force
include("../examples/dynamicStructural/tipSineLoadedCantilever.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicStructural/tipSineLoadedCantilever"))
writedlm("test/newTestDataGenerators/dynamicStructural/tipSineLoadedCantilever/u3_tip.txt", test_data_matrix(u3_tip))
writedlm("test/newTestDataGenerators/dynamicStructural/tipSineLoadedCantilever/F3_root.txt", test_data_matrix(F3_root))
writedlm("test/newTestDataGenerators/dynamicStructural/tipSineLoadedCantilever/M2_root.txt", test_data_matrix(M2_root))

# Dynamic analysis of a cantilever with tip sinusoidal force, with different time integration methods
include("../examples/dynamicStructural/tipSineLoadedCantileverIntegrationMethods.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicStructural/tipSineLoadedCantileverIntegrationMethods"))
writedlm("test/newTestDataGenerators/dynamicStructural/tipSineLoadedCantileverIntegrationMethods/u3_tip.txt", test_data_matrix(u3_tip))
writedlm("test/newTestDataGenerators/dynamicStructural/tipSineLoadedCantileverIntegrationMethods/M2_root.txt", test_data_matrix(M2_root))

# Dynamic analysis of a cantilever with tip sinusoidal force and a tension-only cable
include("../examples/dynamicStructural/tipSineLoadedCantileverWithCable.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/dynamicStructural/tipSineLoadedCantileverWithCable"))
writedlm("test/newTestDataGenerators/dynamicStructural/tipSineLoadedCantileverWithCable/u3_tip.txt", test_data_matrix(u3_tipNorm))
writedlm("test/newTestDataGenerators/dynamicStructural/tipSineLoadedCantileverWithCable/F3_root.txt", test_data_matrix(F3_rootNorm))
writedlm("test/newTestDataGenerators/dynamicStructural/tipSineLoadedCantileverWithCable/M2_root.txt", test_data_matrix(M2_rootNorm))
