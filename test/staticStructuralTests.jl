@testset "Static analysis of an arch under a dead pressure load" begin
    include_example("examples/staticStructural/archUnderDeadPressure.jl")
    @test matches_test_data(mid_u3, "staticStructural/archUnderDeadPressure", "mid_u3.txt")
end

@testset "Static analysis of an arch under a follower pressure load" begin
    include_example("examples/staticStructural/archUnderFollowerPressure.jl")
    @test matches_test_data(mid_u3, "staticStructural/archUnderFollowerPressure", "mid_u3.txt")
end

@testset "Static analysis of a cantilever beam with an axial spring attached between its middle and tip nodes, subjected to an axial tip force" begin
    include_example("examples/staticStructural/axialDoublyAttachedSpringCantilever.jl")
    @test matches_test_data(u1, "staticStructural/axialDoublyAttachedSpringCantilever", "u1.txt")
    @test matches_test_data(F1, "staticStructural/axialDoublyAttachedSpringCantilever", "F1.txt")
end

@testset "Static analysis of a cantilever beam with an axial spring attached between its middle and tip nodes, subjected to an axial tip displacement" begin
    include_example("examples/staticStructural/axialDoublyAttachedSpringCantilever2.jl")
    @test matches_test_data(u1, "staticStructural/axialDoublyAttachedSpringCantilever2", "u1.txt")
    @test matches_test_data(F1, "staticStructural/axialDoublyAttachedSpringCantilever2", "F1.txt")
end

@testset "Static analysis of a biclamped, hinged beam under distributed and concentrated loads" begin
    include_example("examples/staticStructural/biclampedHingedBeam.jl")
    @test matches_test_data(u3Mid, "staticStructural/biclampedHingedBeam", "u3Mid.txt")
    @test matches_test_data(p2Left, "staticStructural/biclampedHingedBeam", "p2Left.txt")
    @test matches_test_data(p2Right, "staticStructural/biclampedHingedBeam", "p2Right.txt")
end

@testset "Static analysis of a cantilever beam bending under self weight" begin
    include_example("examples/staticStructural/cantileverUnderSelfWeight.jl")
    @test matches_test_data(u1, "staticStructural/cantileverUnderSelfWeight", "u1.txt")
    @test matches_test_data(u3, "staticStructural/cantileverUnderSelfWeight", "u3.txt")
    @test matches_test_data(F3, "staticStructural/cantileverUnderSelfWeight", "F3.txt")
    @test matches_test_data(M2, "staticStructural/cantileverUnderSelfWeight", "M2.txt")
end

@testset "Static analysis of a cantilever beam with a tip spring in bending" begin
    include_example("examples/staticStructural/cantileverWithTipSpring.jl")
    @test matches_test_data(u3, "staticStructural/cantileverWithTipSpring", "u3.txt")
    @test matches_test_data(F3, "staticStructural/cantileverWithTipSpring", "F3.txt")
    @test matches_test_data(M2, "staticStructural/cantileverWithTipSpring", "M2.txt")
end

@testset "Static analysis of a clamped beam, with a free flared hinge at the middle, under distributed and concentrated loads" begin
    include_example("examples/staticStructural/clampedFlaredHingedBeam.jl")
    @test matches_test_data(u1, "staticStructural/clampedFlaredHingedBeam", "u1.txt")
    @test matches_test_data(u2, "staticStructural/clampedFlaredHingedBeam", "u2.txt")
    @test matches_test_data(u3, "staticStructural/clampedFlaredHingedBeam", "u3.txt")
    @test matches_test_data(p1, "staticStructural/clampedFlaredHingedBeam", "p1.txt")
    @test matches_test_data(p2, "staticStructural/clampedFlaredHingedBeam", "p2.txt")
    @test matches_test_data(p3, "staticStructural/clampedFlaredHingedBeam", "p3.txt")
    @test matches_test_data(F3, "staticStructural/clampedFlaredHingedBeam", "F3.txt")
    @test matches_test_data(M2, "staticStructural/clampedFlaredHingedBeam", "M2.txt")
    @test matches_test_data(ϕHinge, "staticStructural/clampedFlaredHingedBeam", "phiHinge.txt")
end

@testset "Static analysis of a clamped beam, hinged at the middle with imposed hinge angle, under distributed and concentrated loads" begin
    include_example("examples/staticStructural/clampedHingedBeam.jl")
    @test matches_test_data(u1, "staticStructural/clampedHingedBeam", "u1.txt")
    @test matches_test_data(u3, "staticStructural/clampedHingedBeam", "u3.txt")
    @test matches_test_data(p2, "staticStructural/clampedHingedBeam", "p2.txt")
    @test matches_test_data(F3, "staticStructural/clampedHingedBeam", "F3.txt")
    @test matches_test_data(M2, "staticStructural/clampedHingedBeam", "M2.txt")
end

@testset "Static analysis of a clamped beam rotated in 3D space, hinged at the middle with imposed hinge angle, under distributed and concentrated loads" begin
    include_example("examples/staticStructural/clampedHingedBeamRotated.jl")
    @test matches_test_data(u1, "staticStructural/clampedHingedBeamRotated", "u1.txt")
    @test matches_test_data(u2, "staticStructural/clampedHingedBeamRotated", "u2.txt")
    @test matches_test_data(u3, "staticStructural/clampedHingedBeamRotated", "u3.txt")
    @test matches_test_data(p2_b, "staticStructural/clampedHingedBeamRotated", "p2_b.txt")
    @test matches_test_data(F3, "staticStructural/clampedHingedBeamRotated", "F3.txt")
    @test matches_test_data(M2, "staticStructural/clampedHingedBeamRotated", "M2.txt")
end

@testset "Static analysis of a clamped beam, hinged at the middle, under distributed loads, with varying spring stiffnesses around the hinge" begin
    include_example("examples/staticStructural/clampedHingedBeamSpringRange.jl")
    @test matches_test_data(u1, "staticStructural/clampedHingedBeamSpringRange", "u1.txt")
    @test matches_test_data(u3, "staticStructural/clampedHingedBeamSpringRange", "u3.txt")
    @test matches_test_data(p2, "staticStructural/clampedHingedBeamSpringRange", "p2.txt")
    @test matches_test_data(F3, "staticStructural/clampedHingedBeamSpringRange", "F3.txt")
    @test matches_test_data(M2, "staticStructural/clampedHingedBeamSpringRange", "M2.txt")
end

@testset "Static analysis of a clamped beam with a coasting flared hinge at the middle, under distributed loads" begin
    include_example("examples/staticStructural/coastingFoldingWingtip.jl")
    @test matches_test_data(pHinge, "staticStructural/coastingFoldingWingtip", "pHinge.txt")
end

@testset "Static analysis of composite laminates subjected to tip loads" begin
    include_example("examples/staticStructural/compositeCantileverMD.jl")
    @test matches_test_data(u1_500mm, "staticStructural/compositeCantileverMD", "u1_500mm.txt")
    @test matches_test_data(u2_500mm, "staticStructural/compositeCantileverMD", "u2_500mm.txt")
    @test matches_test_data(u3_500mm, "staticStructural/compositeCantileverMD", "u3_500mm.txt")
end

@testset "Static analysis of a curved cantilever subjected to a tip dead force" begin
    include_example("examples/staticStructural/curvedCantileverDeadLoad.jl")
    @test matches_test_data(tip_u1, "staticStructural/curvedCantileverDeadLoad", "tip_u1.txt")
    @test matches_test_data(tip_u2, "staticStructural/curvedCantileverDeadLoad", "tip_u2.txt")
    @test matches_test_data(tip_u3, "staticStructural/curvedCantileverDeadLoad", "tip_u3.txt")
end

@testset "Static analysis of a curved cantilever subjected to a tip follower force" begin
    include_example("examples/staticStructural/curvedCantileverStaticFollower.jl")
    @test matches_test_data(tip_u1, "staticStructural/curvedCantileverStaticFollower", "tip_u1.txt")
    @test matches_test_data(tip_u2, "staticStructural/curvedCantileverStaticFollower", "tip_u2.txt")
    @test matches_test_data(tip_u3, "staticStructural/curvedCantileverStaticFollower", "tip_u3.txt")
end

@testset "Static analysis of a cantilever with distributed follower force" begin
    include_example("examples/staticStructural/distributedLoadCantilever.jl")
    @test matches_test_data(tip_u1, "staticStructural/distributedLoadCantilever", "tip_u1.txt")
    @test matches_test_data(tip_u3, "staticStructural/distributedLoadCantilever", "tip_u3.txt")
    @test matches_test_data(tip_angle, "staticStructural/distributedLoadCantilever", "tip_angle.txt")
end

@testset "Static analysis of a clamped beam with a driven flared hinge at the middle, under distributed loads" begin
    include_example("examples/staticStructural/drivenFoldingWingtip.jl")
    @test matches_test_data(pHinge, "staticStructural/drivenFoldingWingtip", "pHinge.txt")
end

@testset "Static analysis of a clamped beam with a driven, springed, flared hinge at the middle" begin
    include_example("examples/staticStructural/drivenSpringedHingedBeam.jl")
    @test matches_test_data(pHinge, "staticStructural/drivenSpringedHingedBeam", "pHinge.txt")
end

@testset "Static analysis of a hinged beam subjected to a distributed load" begin
    include_example("examples/staticStructural/hingedBeam.jl")
    @test matches_test_data(u3, "staticStructural/hingedBeam", "u3.txt")
    @test matches_test_data(p2, "staticStructural/hingedBeam", "p2.txt")
    @test matches_test_data(F3, "staticStructural/hingedBeam", "F3.txt")
    @test matches_test_data(M2, "staticStructural/hingedBeam", "M2.txt")
end

@testset "Static analysis of a hinged beam subjected to a distributed load and a rotational spring" begin
    include_example("examples/staticStructural/hingedSpringedBeam.jl")
    @test matches_test_data(u3, "staticStructural/hingedSpringedBeam", "u3.txt")
    @test matches_test_data(p2, "staticStructural/hingedSpringedBeam", "p2.txt")
    @test matches_test_data(F3, "staticStructural/hingedSpringedBeam", "F3.txt")
    @test matches_test_data(M2, "staticStructural/hingedSpringedBeam", "M2.txt")
end

@testset "Static analysis of a T-frame hinged at the connection subjected to a distributed load" begin
    include_example("examples/staticStructural/hingedTFrame.jl")
    @test matches_test_data(u1_beam1, "staticStructural/hingedTFrame", "u1_beam1.txt")
    @test matches_test_data(u1_beam2, "staticStructural/hingedTFrame", "u1_beam2.txt")
    @test matches_test_data(u3_beam1, "staticStructural/hingedTFrame", "u3_beam1.txt")
    @test matches_test_data(u3_beam2, "staticStructural/hingedTFrame", "u3_beam2.txt")
    @test matches_test_data(F3_beam1, "staticStructural/hingedTFrame", "F3_beam1.txt")
    @test matches_test_data(M2_beam1, "staticStructural/hingedTFrame", "M2_beam1.txt")
end

@testset "Static analysis of the Lee frame (a right-angled frame) with a dead load" begin
    include_example("examples/staticStructural/LeeFrameDeadLoad.jl")
    @test matches_test_data(u1_atForce, "staticStructural/LeeFrameDeadLoad", "u1_atForce.txt")
    @test matches_test_data(u3_atForce, "staticStructural/LeeFrameDeadLoad", "u3_atForce.txt")
end

@testset "Static analysis of the Lee frame (a right-angled frame) with a follower load" begin
    include_example("examples/staticStructural/LeeFrameFollowerLoad.jl")
    @test matches_test_data(u1_atForce, "staticStructural/LeeFrameFollowerLoad", "u1_atForce.txt")
    @test matches_test_data(u3_atForce, "staticStructural/LeeFrameFollowerLoad", "u3_atForce.txt")
end

@testset "Static analysis of two cantilevers connected at their tips by nonlinear doubly-attached springs in several directions, with and without rigid-body rotations of the assembly" begin
    include_example("examples/staticStructural/nonlinearDoublyAttachedSprings.jl")
    @test matches_test_data(tipDisp, "staticStructural/nonlinearDoublyAttachedSprings", "tipDisp.txt")
    @test matches_test_data(clampLoads, "staticStructural/nonlinearDoublyAttachedSprings", "clampLoads.txt")
end

@testset "Static analysis of the pure bending test of the Pazy wing" begin
    include_example("examples/staticStructural/PazyWingBendingTest.jl")
    @test matches_test_data(tip_OOP, "staticStructural/PazyWingBendingTest", "tip_OOP.txt")
end

@testset "Static analysis of the coupled torsion-bending test of the Pazy wing" begin
    include_example("examples/staticStructural/PazyWingTorsionTest.jl")
    @test matches_test_data(tip_OOP, "staticStructural/PazyWingTorsionTest", "tip_OOP.txt")
    @test matches_test_data(tip_twist, "staticStructural/PazyWingTorsionTest", "tip_twist.txt")
end

@testset "Static analysis of a mid-loaded arch pinned at one end and clamped at the other" begin
    include_example("examples/staticStructural/pinnedClampedArch.jl")
    @test matches_test_data(u1_atForce, "staticStructural/pinnedClampedArch", "u1.txt")
    @test matches_test_data(u3_atForce, "staticStructural/pinnedClampedArch", "u3.txt")
end

@testset "Static analysis of a right-angled frame under a tip transverse follower force" begin
    include_example("examples/staticStructural/rightAngledFrame.jl")
    @test matches_test_data(tip_u1, "staticStructural/rightAngledFrame", "tip_u1.txt")
    @test matches_test_data(tip_u3, "staticStructural/rightAngledFrame", "tip_u3.txt")
    @test matches_test_data(tip_angle, "staticStructural/rightAngledFrame", "tip_angle.txt")
end

@testset "Static analysis of a right-angled frame under a 'buckling' load" begin
    include_example("examples/staticStructural/rightAngledFrameBuckling.jl")
    @test matches_test_data(tip_u2, "staticStructural/rightAngledFrameBuckling", "tip_u2.txt")
end

@testset "Static analysis of a L-frame with a doubly-attached spring and tip load" begin
    include_example("examples/staticStructural/springedLFrame.jl")
    @test matches_test_data(u3_b, "staticStructural/springedLFrame", "u3_b.txt")
end

@testset "Static analysis of a swept-back clamped beam with a driven flared folding wingtip" begin
    include_example("examples/staticStructural/sweptBackFFWTWing.jl")
    @test matches_test_data(ϕHinge, "staticStructural/sweptBackFFWTWing", "phi.txt")
    @test matches_test_data(hingeMoment, "staticStructural/sweptBackFFWTWing", "hingeMoment.txt")
end

@testset "Static analysis of a semi-circular arch with tangential follower force" begin
    include_example("examples/staticStructural/tangentiallyForcedArch.jl")
    @test matches_test_data(tip_u1, "staticStructural/tangentiallyForcedArch", "tip_u1.txt")
    @test matches_test_data(tip_u3, "staticStructural/tangentiallyForcedArch", "tip_u3.txt")
    @test matches_test_data(tip_angle, "staticStructural/tangentiallyForcedArch", "tip_angle.txt")
end

@testset "Static analysis cantilever with tip follower transverse force" begin
    include_example("examples/staticStructural/tipFollowerForceCantilever.jl")
    @test matches_test_data(tip_u1, "staticStructural/tipFollowerForceCantilever", "tip_u1.txt")
    @test matches_test_data(tip_u3, "staticStructural/tipFollowerForceCantilever", "tip_u3.txt")
end

@testset "Static analysis of a cantilever with tip follower transverse force (force split over 2 BCs)" begin
    include_example("examples/staticStructural/tipFollowerForceCantilever2.jl")
    @test matches_test_data(tip_u1, "staticStructural/tipFollowerForceCantilever2", "tip_u1.txt")
    @test matches_test_data(tip_u3, "staticStructural/tipFollowerForceCantilever2", "tip_u3.txt")
end

@testset "Static analysis of a cantilever with tip moment" begin
    include_example("examples/staticStructural/tipMomentCantilever.jl")
    @test matches_test_data(tip_u1, "staticStructural/tipMomentCantilever", "tip_u1.txt")
    @test matches_test_data(tip_u3, "staticStructural/tipMomentCantilever", "tip_u3.txt")
    @test matches_test_data(tip_angle, "staticStructural/tipMomentCantilever", "tip_angle.txt")
end

@testset "Static analysis of two cantilevers connected at their tips by a transverse doubly-attached spring, at several locations of the spring, against the analytical solution" begin
    include_example("examples/staticStructural/transverseDoublyAttachedSpringCantilevers.jl")
    @test matches_test_data(springDisp, "staticStructural/transverseDoublyAttachedSpringCantilevers", "springDisp.txt")
    @test matches_test_data(tipsDisp, "staticStructural/transverseDoublyAttachedSpringCantilevers", "tipDisp.txt")
    @test matches_test_data(tipsTwist, "staticStructural/transverseDoublyAttachedSpringCantilevers", "tipTwist.txt")
    @test ϵAnalytical < 1e-6
end

@testset "Static analysis of a semi-circular arch with transverse follower force" begin
    include_example("examples/staticStructural/transverselyForcedArch.jl")
    @test matches_test_data(tip_u1, "staticStructural/transverselyForcedArch", "tip_u1.txt")
    @test matches_test_data(tip_u3, "staticStructural/transverselyForcedArch", "tip_u3.txt")
    @test matches_test_data(tip_angle, "staticStructural/transverselyForcedArch", "tip_angle.txt")
end

@testset "Static analysis of a beam with triangular distributed load" begin
    include_example("examples/staticStructural/triangleLoadBeam.jl")
    @test matches_test_data(u3, "staticStructural/triangleLoadBeam", "u3.txt")
    @test matches_test_data(F3, "staticStructural/triangleLoadBeam", "F3.txt")
    @test matches_test_data(M2, "staticStructural/triangleLoadBeam", "M2.txt")
end

@testset "Static analysis of a cantilever beam with a torsional spring attached between its middle and tip nodes, subjected to a torsional tip moment" begin
    include_example("examples/staticStructural/twistDoublyAttachedSpringCantilever.jl")
    @test matches_test_data(p1, "staticStructural/twistDoublyAttachedSpringCantilever", "p1.txt")
    @test matches_test_data(M1, "staticStructural/twistDoublyAttachedSpringCantilever", "M1.txt")
end

@testset "Static analysis of a beam with a central hinge driven by opposite end moments" begin
    include_example("examples/staticStructural/drivenSpringedHingedBeam2.jl")
    # Self-comparison
    @test matches_test_data(u3, "staticStructural/drivenSpringedHingedBeam2", "u3.txt")
    @test matches_test_data(p2, "staticStructural/drivenSpringedHingedBeam2", "p2.txt")
    @test matches_test_data(M2, "staticStructural/drivenSpringedHingedBeam2", "M2.txt")
end

@testset "Static analysis of a swept cantilever under tip transverse loading" begin
    include_example("examples/staticStructural/sweptBeamStaticLoading.jl")
    # Self-comparison
    @test matches_test_data(tip_u3, "staticStructural/sweptBeamStaticLoading", "tip_u3.txt")
    @test matches_test_data(tip_twist, "staticStructural/sweptBeamStaticLoading", "tip_twist.txt")
end