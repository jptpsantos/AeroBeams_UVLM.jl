@testset "Dynamic analysis of the axial vibration of a beam under a traction force applied suddenly" begin
    include_example("examples/dynamicStructural/axialTractionCantilever.jl")
    @test matches_test_data(u1_08, "dynamicStructural/axialTractionCantilever", "u1_08.txt")
    @test matches_test_data(u1_10, "dynamicStructural/axialTractionCantilever", "u1_10.txt")
end

@testset "Dynamic analysis of the free response of a beam clamped at both ends and subjected to an initial displacement profile" begin
    include_example("examples/dynamicStructural/biclampedBeam.jl")
    @test matches_test_data(u3_mid, "dynamicStructural/biclampedBeam", "u3_mid.txt")
    @test matches_test_data(V3_mid, "dynamicStructural/biclampedBeam", "V3_mid.txt")
    @test matches_test_data(Vdot3_mid, "dynamicStructural/biclampedBeam", "Vdot3_mid.txt")
    @test matches_test_data(θ2_quarter, "dynamicStructural/biclampedBeam", "theta2_quarter.txt")
    @test matches_test_data(Ω2_quarter, "dynamicStructural/biclampedBeam", "Omega2_quarter.txt")
    @test matches_test_data(Ωdot2_quarter, "dynamicStructural/biclampedBeam", "Omegadot2_quarter.txt")
end

@testset "Dynamic analysis of the free decay of a cantilever with stiffness-proportional structural damping" begin
    include_example("examples/dynamicStructural/cantileverStructuralDampingFreeDecay.jl")
    @test abs(ratio / ratioExpected - 1) < 1e-3
    @test matches_test_data(u3Tip, "dynamicStructural/cantileverStructuralDampingFreeDecay", "u3_tip.txt")
end

@testset "Dynamic analysis of a slanted double pendulum released from rest" begin
    include_example("examples/dynamicStructural/ChengSlantedDoublePendulum.jl")
    @test matches_test_data(u1_tip, "dynamicStructural/ChengSlantedDoublePendulum", "u1_tip.txt")
    @test matches_test_data(u2_tip, "dynamicStructural/ChengSlantedDoublePendulum", "u2_tip.txt")
    @test matches_test_data(u3_tip, "dynamicStructural/ChengSlantedDoublePendulum", "u3_tip.txt")
end

@testset "Dynamic analysis of a composite cantilever beam under a tip sinusoidal load" begin
    include_example("examples/dynamicStructural/compositeCantilever.jl")
    @test matches_test_data(u1_tip, "dynamicStructural/compositeCantilever", "u1_tip.txt")
    @test matches_test_data(u2_tip, "dynamicStructural/compositeCantilever", "u2_tip.txt")
    @test matches_test_data(u3_tip, "dynamicStructural/compositeCantilever", "u3_tip.txt")
    @test matches_test_data(p1_tip, "dynamicStructural/compositeCantilever", "p1_tip.txt")
    @test matches_test_data(p2_tip, "dynamicStructural/compositeCantilever", "p2_tip.txt")
    @test matches_test_data(p3_tip, "dynamicStructural/compositeCantilever", "p3_tip.txt")
    @test matches_test_data(F1_root, "dynamicStructural/compositeCantilever", "F1_root.txt")
    @test matches_test_data(F2_root, "dynamicStructural/compositeCantilever", "F2_root.txt")
    @test matches_test_data(F3_root, "dynamicStructural/compositeCantilever", "F3_root.txt")
    @test matches_test_data(M1_root, "dynamicStructural/compositeCantilever", "M1_root.txt")
    @test matches_test_data(M2_root, "dynamicStructural/compositeCantilever", "M2_root.txt")
    @test matches_test_data(M3_root, "dynamicStructural/compositeCantilever", "M3_root.txt")
end

@testset "Dynamic analysis of a curved cantilever subjected to a tip follower force" begin
    include_example("examples/dynamicStructural/curvedCantileverDynamicFollower.jl")
    @test matches_test_data(u1_tip, "dynamicStructural/curvedCantileverDynamicFollower", "u1_tip.txt")
    @test matches_test_data(u2_tip, "dynamicStructural/curvedCantileverDynamicFollower", "u2_tip.txt")
    @test matches_test_data(u3_tip, "dynamicStructural/curvedCantileverDynamicFollower", "u3_tip.txt")
    @test matches_test_data(F1_root, "dynamicStructural/curvedCantileverDynamicFollower", "F1_root.txt")
    @test matches_test_data(F2_root, "dynamicStructural/curvedCantileverDynamicFollower", "F2_root.txt")
    @test matches_test_data(F3_root, "dynamicStructural/curvedCantileverDynamicFollower", "F3_root.txt")
    @test matches_test_data(M1_root, "dynamicStructural/curvedCantileverDynamicFollower", "M1_root.txt")
    @test matches_test_data(M2_root, "dynamicStructural/curvedCantileverDynamicFollower", "M2_root.txt")
    @test matches_test_data(M3_root, "dynamicStructural/curvedCantileverDynamicFollower", "M3_root.txt")
end

@testset "Dynamic analysis of a double pendulum released from rest" begin
    include_example("examples/dynamicStructural/doublePendulum.jl")
    @test matches_test_data(u1_hinge, "dynamicStructural/doublePendulum", "u1_hinge.txt")
    @test matches_test_data(u3_hinge, "dynamicStructural/doublePendulum", "u3_hinge.txt")
    @test matches_test_data(u1_tip, "dynamicStructural/doublePendulum", "u1_tip.txt")
    @test matches_test_data(u3_tip, "dynamicStructural/doublePendulum", "u3_tip.txt")
end

@testset "Dynamic analysis of a right-angled frame subjected to an out-of-plane force" begin
    include_example("examples/dynamicStructural/elbowFrame.jl")
    @test matches_test_data(u3_elbow, "dynamicStructural/elbowFrame", "u3_elbow.txt")
    @test matches_test_data(u3_tip, "dynamicStructural/elbowFrame", "u3_tip.txt")
end

@testset "Dynamic analysis of a flexible beam subjected to loads yielding two-dimensional motion" begin
    include_example("examples/dynamicStructural/flyingFlexibleBeam2D.jl")
    @test matches_test_data(u1_tip, "dynamicStructural/flyingFlexibleBeam2D", "u1_tip.txt")
    @test matches_test_data(u3_tip, "dynamicStructural/flyingFlexibleBeam2D", "u3_tip.txt")
end

@testset "Dynamic analysis of a hinged beam in free flight" begin
    include_example("examples/dynamicStructural/flyingScissors.jl")
    @test matches_test_data(u1_tipA, "dynamicStructural/flyingScissors", "u1_tipA.txt")
    @test matches_test_data(u3_tipA, "dynamicStructural/flyingScissors", "u3_tipA.txt")
    @test matches_test_data(u1_tipB, "dynamicStructural/flyingScissors", "u1_tipB.txt")
    @test matches_test_data(u3_tipB, "dynamicStructural/flyingScissors", "u3_tipB.txt")
    @test matches_test_data(u1_hinge, "dynamicStructural/flyingScissors", "u1_hinge.txt")
    @test matches_test_data(u3_hinge, "dynamicStructural/flyingScissors", "u3_hinge.txt")
end

@testset "Dynamic analysis of a very flexible beam subjected to loads yielding two-dimensional motion" begin
    include_example("examples/dynamicStructural/flyingSpaghetti2D.jl")
    @test matches_test_data(u1_tip, "dynamicStructural/flyingSpaghetti2D", "u1_tip.txt")
    @test matches_test_data(u3_tip, "dynamicStructural/flyingSpaghetti2D", "u3_tip.txt")
end

@testset "Dynamic analysis of a very flexible beam subjected to loads yielding tri-dimensional motion" begin
    include_example("examples/dynamicStructural/flyingSpaghetti3D.jl")
    @test matches_test_data(u1_tip, "dynamicStructural/flyingSpaghetti3D", "u1_tip.txt")
    @test matches_test_data(u2_tip, "dynamicStructural/flyingSpaghetti3D", "u2_tip.txt")
    @test matches_test_data(u3_tip, "dynamicStructural/flyingSpaghetti3D", "u3_tip.txt")
end

@testset "Dynamic analysis of the free response of a beam subjected to initial displacement and velocity profiles" begin
    include_example("examples/dynamicStructural/initialDispAndVelBeam.jl")
    @test matches_test_data(u3_quarter, "dynamicStructural/initialDispAndVelBeam", "u3_quarter.txt")
    @test matches_test_data(V3_quarter, "dynamicStructural/initialDispAndVelBeam", "V3_quarter.txt")
    @test matches_test_data(Vdot3_quarter, "dynamicStructural/initialDispAndVelBeam", "Vdot3_quarter.txt")
    @test matches_test_data(θ2_root, "dynamicStructural/initialDispAndVelBeam", "theta2_root.txt")
    @test matches_test_data(Ω2_mid, "dynamicStructural/initialDispAndVelBeam", "Omega2_mid.txt")
    @test matches_test_data(Ωdot2_mid, "dynamicStructural/initialDispAndVelBeam", "Omegadot2_mid.txt")
end

@testset "Dynamic analysis of the free response of a beam subjected to an initial displacement profile" begin
    include_example("examples/dynamicStructural/initialDisplacementBeam.jl")
    @test matches_test_data(u3_quarter, "dynamicStructural/initialDisplacementBeam", "u3_quarter.txt")
    @test matches_test_data(V3_quarter, "dynamicStructural/initialDisplacementBeam", "V3_quarter.txt")
    @test matches_test_data(Vdot3_quarter, "dynamicStructural/initialDisplacementBeam", "Vdot3_quarter.txt")
    @test matches_test_data(θ2_root, "dynamicStructural/initialDisplacementBeam", "theta2_root.txt")
    @test matches_test_data(Ω2_mid, "dynamicStructural/initialDisplacementBeam", "Omega2_mid.txt")
    @test matches_test_data(Ωdot2_mid, "dynamicStructural/initialDisplacementBeam", "Omegadot2_mid.txt")
end

@testset "Dynamic analysis of the free response of a beam subjected to an initial velocity profile" begin
    include_example("examples/dynamicStructural/initialVelocityBeam.jl")
    @test matches_test_data(u3_quarter, "dynamicStructural/initialVelocityBeam", "u3_quarter.txt")
    @test matches_test_data(V3_quarter, "dynamicStructural/initialVelocityBeam", "V3_quarter.txt")
    @test matches_test_data(Vdot3_quarter, "dynamicStructural/initialVelocityBeam", "Vdot3_quarter.txt")
    @test matches_test_data(θ2_root, "dynamicStructural/initialVelocityBeam", "theta2_root.txt")
    @test matches_test_data(Ω2_mid, "dynamicStructural/initialVelocityBeam", "Omega2_mid.txt")
    @test matches_test_data(Ωdot2_mid, "dynamicStructural/initialVelocityBeam", "Omegadot2_mid.txt")
end

@testset "Dynamic analysis of joined beams under load" begin
    include_example("examples/dynamicStructural/joinedBeams.jl")
    @test matches_test_data(u1, "dynamicStructural/joinedBeams", "u1.txt")
    @test matches_test_data(u2, "dynamicStructural/joinedBeams", "u2.txt")
    @test matches_test_data(u3, "dynamicStructural/joinedBeams", "u3.txt")
end

@testset "Dynamic analysis of a pinned robot arm driven by a couple moment" begin
    include_example("examples/dynamicStructural/momentDrivenRobotArm.jl")
    @test matches_test_data(u1_tip, "dynamicStructural/momentDrivenRobotArm", "u1_tip.txt")
    @test matches_test_data(u3_tip, "dynamicStructural/momentDrivenRobotArm", "u3_tip.txt")
end

@testset "Dynamic analysis of a pendulum released from rest" begin
    include_example("examples/dynamicStructural/pendulum.jl")
    @test matches_test_data(u1_tip, "dynamicStructural/pendulum", "u1_tip.txt")
    @test matches_test_data(u3_tip, "dynamicStructural/pendulum", "u3_tip.txt")
end

@testset "Dynamic analysis of a clamped beam with root sinusoidal oscillation (ω = 9 Hz ≈ 2nd bending mode)" begin
    include_example("examples/dynamicStructural/rootExcitationBeam.jl")
    @test matches_test_data(u3b_root, "dynamicStructural/rootExcitationBeam", "u3b_root.txt")
    @test matches_test_data(u3b_tip, "dynamicStructural/rootExcitationBeam", "u3b_tip.txt")
    @test matches_test_data(V3_root, "dynamicStructural/rootExcitationBeam", "V3_root.txt")
    @test matches_test_data(V3_tip, "dynamicStructural/rootExcitationBeam", "V3_tip.txt")
end

@testset "Dynamic analysis of a rotary shaft with specified rotation" begin
    include_example("examples/dynamicStructural/rotaryShaft.jl")
    @test matches_test_data(pNum, "dynamicStructural/rotaryShaft", "pNum.txt")
    @test matches_test_data(pdotNum, "dynamicStructural/rotaryShaft", "pdotNum.txt")
    @test matches_test_data(ΩNum, "dynamicStructural/rotaryShaft", "OmegaNum.txt")
    @test matches_test_data(MNum, "dynamicStructural/rotaryShaft", "MNum.txt")
end

@testset "Dynamic analysis of an articulated robot arm driven by specified rotation" begin
    include_example("examples/dynamicStructural/rotationDrivenArticulatedRobotArm.jl")
    @test matches_test_data(u1_tip, "dynamicStructural/rotationDrivenArticulatedRobotArm", "u1_tip.txt")
    @test matches_test_data(u3_tip, "dynamicStructural/rotationDrivenArticulatedRobotArm", "u3_tip.txt")
    @test matches_test_data(u1_hinge, "dynamicStructural/rotationDrivenArticulatedRobotArm", "u1_hinge.txt")
    @test matches_test_data(u3_hinge, "dynamicStructural/rotationDrivenArticulatedRobotArm", "u3_hinge.txt")
end

@testset "Dynamic analysis of a robot arm driven by specified rotation" begin
    include_example("examples/dynamicStructural/rotationDrivenRobotArm.jl")
    @test matches_test_data(u1_tip, "dynamicStructural/rotationDrivenRobotArm", "u1_tip.txt")
    @test matches_test_data(u3_tip, "dynamicStructural/rotationDrivenRobotArm", "u3_tip.txt")
end

@testset "Dynamic analysis of the spin-up maneuver of a robot arm" begin
    include_example("examples/dynamicStructural/spinupRobotArm.jl")
    @test matches_test_data(u1_tip, "dynamicStructural/spinupRobotArm", "u1_tip.txt")
    @test matches_test_data(u2_tip, "dynamicStructural/spinupRobotArm", "u2_tip.txt")
    @test matches_test_data(θ3_root, "dynamicStructural/spinupRobotArm", "theta3_root.txt")
end

@testset "Dynamic analysis of a cantilever with linearly increasing tip moment" begin
    include_example("examples/dynamicStructural/tipMomentCantileverDynamic.jl")
    @test matches_test_data(κ2, "dynamicStructural/tipMomentCantileverDynamic", "kappa2.txt")
end

@testset "Dynamic analysis of a pendulum with tip mass" begin
    include_example("examples/dynamicStructural/tipPendulum.jl")
    @test matches_test_data(u1_tip, "dynamicStructural/tipPendulum", "u1_tip.txt")
    @test matches_test_data(u3_tip, "dynamicStructural/tipPendulum", "u3_tip.txt")
end

@testset "Dynamic analysis of a cantilever with tip sinusoidal force" begin
    include_example("examples/dynamicStructural/tipSineLoadedCantilever.jl")
    @test matches_test_data(u3_tip, "dynamicStructural/tipSineLoadedCantilever", "u3_tip.txt")
    @test matches_test_data(F3_root, "dynamicStructural/tipSineLoadedCantilever", "F3_root.txt")
    @test matches_test_data(M2_root, "dynamicStructural/tipSineLoadedCantilever", "M2_root.txt")
end

@testset "Dynamic analysis of a cantilever with tip sinusoidal force, with different time integration methods" begin
    include_example("examples/dynamicStructural/tipSineLoadedCantileverIntegrationMethods.jl")
    @test matches_test_data(u3_tip, "dynamicStructural/tipSineLoadedCantileverIntegrationMethods", "u3_tip.txt")
    @test matches_test_data(M2_root, "dynamicStructural/tipSineLoadedCantileverIntegrationMethods", "M2_root.txt")
end

@testset "Dynamic analysis of a cantilever with tip sinusoidal force and a tension-only cable" begin
    include_example("examples/dynamicStructural/tipSineLoadedCantileverWithCable.jl")
    @test matches_test_data(u3_tipNorm, "dynamicStructural/tipSineLoadedCantileverWithCable", "u3_tipNorm.txt")
    @test matches_test_data(F3_rootNorm, "dynamicStructural/tipSineLoadedCantileverWithCable", "F3_rootNorm.txt")
    @test matches_test_data(M2_rootNorm, "dynamicStructural/tipSineLoadedCantileverWithCable", "M2_rootNorm.txt")
end