# Static analysis of an arch under a dead pressure load
include("../examples/staticStructural/archUnderDeadPressure.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/archUnderDeadPressure"))
writedlm("test/newTestDataGenerators/staticStructural/archUnderDeadPressure/mid_u3.txt", test_data_matrix(mid_u3))

# Static analysis of an arch under a follower pressure load
include("../examples/staticStructural/archUnderFollowerPressure.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/archUnderFollowerPressure"))
writedlm("test/newTestDataGenerators/staticStructural/archUnderFollowerPressure/mid_u3.txt", test_data_matrix(mid_u3))

# Static analysis of a cantilever beam with an axial spring attached between its middle and tip nodes, subjected to an axial tip force
include("../examples/staticStructural/axialDoublyAttachedSpringCantilever.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/axialDoublyAttachedSpringCantilever"))
writedlm("test/newTestDataGenerators/staticStructural/axialDoublyAttachedSpringCantilever/u1.txt", test_data_matrix(u1))
writedlm("test/newTestDataGenerators/staticStructural/axialDoublyAttachedSpringCantilever/F1.txt", test_data_matrix(F1))

# Static analysis of a cantilever beam with an axial spring attached between its middle and tip nodes, subjected to an axial tip displacement
include("../examples/staticStructural/axialDoublyAttachedSpringCantilever2.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/axialDoublyAttachedSpringCantilever2"))
writedlm("test/newTestDataGenerators/staticStructural/axialDoublyAttachedSpringCantilever2/u1.txt", test_data_matrix(u1))
writedlm("test/newTestDataGenerators/staticStructural/axialDoublyAttachedSpringCantilever2/F1.txt", test_data_matrix(F1))

# Static analysis of a biclamped, hinged beam under distributed and concentrated loads
include("../examples/staticStructural/biclampedHingedBeam.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/biclampedHingedBeam"))
writedlm("test/newTestDataGenerators/staticStructural/biclampedHingedBeam/u3Mid.txt", test_data_matrix(u3Mid))
writedlm("test/newTestDataGenerators/staticStructural/biclampedHingedBeam/p2Left.txt", test_data_matrix(p2Left))
writedlm("test/newTestDataGenerators/staticStructural/biclampedHingedBeam/p2Right.txt", test_data_matrix(p2Right))

# Static analysis of a cantilever beam bending under self weight
include("../examples/staticStructural/cantileverUnderSelfWeight.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/cantileverUnderSelfWeight"))
writedlm("test/newTestDataGenerators/staticStructural/cantileverUnderSelfWeight/u1.txt", test_data_matrix(u1))
writedlm("test/newTestDataGenerators/staticStructural/cantileverUnderSelfWeight/u3.txt", test_data_matrix(u3))
writedlm("test/newTestDataGenerators/staticStructural/cantileverUnderSelfWeight/F3.txt", test_data_matrix(F3))
writedlm("test/newTestDataGenerators/staticStructural/cantileverUnderSelfWeight/M2.txt", test_data_matrix(M2))

# Static analysis of a cantilever beam with a tip spring in bending
include("../examples/staticStructural/cantileverWithTipSpring.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/cantileverWithTipSpring"))
writedlm("test/newTestDataGenerators/staticStructural/cantileverWithTipSpring/u3.txt", test_data_matrix(u3))
writedlm("test/newTestDataGenerators/staticStructural/cantileverWithTipSpring/F3.txt", test_data_matrix(F3))
writedlm("test/newTestDataGenerators/staticStructural/cantileverWithTipSpring/M2.txt", test_data_matrix(M2))

# Static analysis of a clamped beam, with a free flared hinge at the middle, under distributed and concentrated loads
include("../examples/staticStructural/clampedFlaredHingedBeam.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/clampedFlaredHingedBeam"))
writedlm("test/newTestDataGenerators/staticStructural/clampedFlaredHingedBeam/u1.txt", test_data_matrix(u1))
writedlm("test/newTestDataGenerators/staticStructural/clampedFlaredHingedBeam/u2.txt", test_data_matrix(u2))
writedlm("test/newTestDataGenerators/staticStructural/clampedFlaredHingedBeam/u3.txt", test_data_matrix(u3))
writedlm("test/newTestDataGenerators/staticStructural/clampedFlaredHingedBeam/p1.txt", test_data_matrix(p1))
writedlm("test/newTestDataGenerators/staticStructural/clampedFlaredHingedBeam/p2.txt", test_data_matrix(p2))
writedlm("test/newTestDataGenerators/staticStructural/clampedFlaredHingedBeam/p3.txt", test_data_matrix(p3))
writedlm("test/newTestDataGenerators/staticStructural/clampedFlaredHingedBeam/F3.txt", test_data_matrix(F3))
writedlm("test/newTestDataGenerators/staticStructural/clampedFlaredHingedBeam/M2.txt", test_data_matrix(M2))
writedlm("test/newTestDataGenerators/staticStructural/clampedFlaredHingedBeam/phiHinge.txt", test_data_matrix(ϕHinge))

# Static analysis of a clamped beam, hinged at the middle with imposed hinge angle, under distributed and concentrated loads
include("../examples/staticStructural/clampedHingedBeam.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/clampedHingedBeam"))
writedlm("test/newTestDataGenerators/staticStructural/clampedHingedBeam/u1.txt", test_data_matrix(u1))
writedlm("test/newTestDataGenerators/staticStructural/clampedHingedBeam/u3.txt", test_data_matrix(u3))
writedlm("test/newTestDataGenerators/staticStructural/clampedHingedBeam/p2.txt", test_data_matrix(p2))
writedlm("test/newTestDataGenerators/staticStructural/clampedHingedBeam/F3.txt", test_data_matrix(F3))
writedlm("test/newTestDataGenerators/staticStructural/clampedHingedBeam/M2.txt", test_data_matrix(M2))

# Static analysis of a clamped beam rotated in 3D space, hinged at the middle with imposed hinge angle, under distributed and concentrated loads
include("../examples/staticStructural/clampedHingedBeamRotated.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/clampedHingedBeamRotated"))
writedlm("test/newTestDataGenerators/staticStructural/clampedHingedBeamRotated/u1.txt", test_data_matrix(u1))
writedlm("test/newTestDataGenerators/staticStructural/clampedHingedBeamRotated/u2.txt", test_data_matrix(u2))
writedlm("test/newTestDataGenerators/staticStructural/clampedHingedBeamRotated/u3.txt", test_data_matrix(u3))
writedlm("test/newTestDataGenerators/staticStructural/clampedHingedBeamRotated/p2_b.txt", test_data_matrix(p2_b))
writedlm("test/newTestDataGenerators/staticStructural/clampedHingedBeamRotated/F3.txt", test_data_matrix(F3))
writedlm("test/newTestDataGenerators/staticStructural/clampedHingedBeamRotated/M2.txt", test_data_matrix(M2))

# Static analysis of a clamped beam, hinged at the middle, under distributed loads, with varying spring stiffnesses around the hinge
include("../examples/staticStructural/clampedHingedBeamSpringRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/clampedHingedBeamSpringRange"))
writedlm("test/newTestDataGenerators/staticStructural/clampedHingedBeamSpringRange/u1.txt", test_data_matrix(u1))
writedlm("test/newTestDataGenerators/staticStructural/clampedHingedBeamSpringRange/u3.txt", test_data_matrix(u3))
writedlm("test/newTestDataGenerators/staticStructural/clampedHingedBeamSpringRange/p2.txt", test_data_matrix(p2))
writedlm("test/newTestDataGenerators/staticStructural/clampedHingedBeamSpringRange/F3.txt", test_data_matrix(F3))
writedlm("test/newTestDataGenerators/staticStructural/clampedHingedBeamSpringRange/M2.txt", test_data_matrix(M2))

# Static analysis of a clamped beam with a coasting flared hinge at the middle, under distributed loads
include("../examples/staticStructural/coastingFoldingWingtip.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/coastingFoldingWingtip"))
writedlm("test/newTestDataGenerators/staticStructural/coastingFoldingWingtip/pHinge.txt", test_data_matrix(pHinge))

# Static analysis of a clamped beam with a driven flared hinge at the middle, under distributed loads
include("../examples/staticStructural/drivenFoldingWingtip.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/drivenFoldingWingtip"))
writedlm("test/newTestDataGenerators/staticStructural/drivenFoldingWingtip/pHinge.txt", test_data_matrix(pHinge))

# Static analysis of a clamped beam with a driven, springed, flared hinge at the middle
include("../examples/staticStructural/drivenSpringedHingedBeam.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/drivenSpringedHingedBeam"))
writedlm("test/newTestDataGenerators/staticStructural/drivenSpringedHingedBeam/pHinge.txt", test_data_matrix(pHinge))

# Static analysis of composite laminates subjected to tip loads
include("../examples/staticStructural/compositeCantileverMD.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/compositeCantileverMD"))
writedlm(string("test/newTestDataGenerators/staticStructural/compositeCantileverMD/u1_500mm.txt"), test_data_matrix(u1_500mm))
writedlm(string("test/newTestDataGenerators/staticStructural/compositeCantileverMD/u2_500mm.txt"), test_data_matrix(u2_500mm))
writedlm(string("test/newTestDataGenerators/staticStructural/compositeCantileverMD/u3_500mm.txt"), test_data_matrix(u3_500mm))

# Static analysis of a curved cantilever subjected to a tip dead force
include("../examples/staticStructural/curvedCantileverDeadLoad.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/curvedCantileverDeadLoad"))
writedlm("test/newTestDataGenerators/staticStructural/curvedCantileverDeadLoad/tip_u1.txt", test_data_matrix(tip_u1))
writedlm("test/newTestDataGenerators/staticStructural/curvedCantileverDeadLoad/tip_u2.txt", test_data_matrix(tip_u2))
writedlm("test/newTestDataGenerators/staticStructural/curvedCantileverDeadLoad/tip_u3.txt", test_data_matrix(tip_u3))

# Static analysis of a curved cantilever subjected to a tip follower force
include("../examples/staticStructural/curvedCantileverStaticFollower.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/curvedCantileverStaticFollower"))
writedlm("test/newTestDataGenerators/staticStructural/curvedCantileverStaticFollower/tip_u1.txt", test_data_matrix(tip_u1))
writedlm("test/newTestDataGenerators/staticStructural/curvedCantileverStaticFollower/tip_u2.txt", test_data_matrix(tip_u2))
writedlm("test/newTestDataGenerators/staticStructural/curvedCantileverStaticFollower/tip_u3.txt", test_data_matrix(tip_u3))

# Static analysis of a cantilever with distributed follower force
include("../examples/staticStructural/distributedLoadCantilever.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/distributedLoadCantilever"))
writedlm("test/newTestDataGenerators/staticStructural/distributedLoadCantilever/tip_u1.txt", test_data_matrix(tip_u1))
writedlm("test/newTestDataGenerators/staticStructural/distributedLoadCantilever/tip_u3.txt", test_data_matrix(tip_u3))
writedlm("test/newTestDataGenerators/staticStructural/distributedLoadCantilever/tip_angle.txt", test_data_matrix(tip_angle))

# Static analysis of a hinged beam subjected to a distributed load
include("../examples/staticStructural/hingedBeam.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/hingedBeam"))
writedlm("test/newTestDataGenerators/staticStructural/hingedBeam/u3.txt", test_data_matrix(u3))
writedlm("test/newTestDataGenerators/staticStructural/hingedBeam/p2.txt", test_data_matrix(p2))
writedlm("test/newTestDataGenerators/staticStructural/hingedBeam/F3.txt", test_data_matrix(F3))
writedlm("test/newTestDataGenerators/staticStructural/hingedBeam/M2.txt", test_data_matrix(M2))

# Static analysis of a hinged beam subjected to a distributed load and a rotational spring
include("../examples/staticStructural/hingedSpringedBeam.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/hingedSpringedBeam"))
writedlm("test/newTestDataGenerators/staticStructural/hingedSpringedBeam/u3.txt", test_data_matrix(u3))
writedlm("test/newTestDataGenerators/staticStructural/hingedSpringedBeam/p2.txt", test_data_matrix(p2))
writedlm("test/newTestDataGenerators/staticStructural/hingedSpringedBeam/F3.txt", test_data_matrix(F3))
writedlm("test/newTestDataGenerators/staticStructural/hingedSpringedBeam/M2.txt", test_data_matrix(M2))

# Static analysis of a T-frame hinged at the connection subjected to a distributed load
include("../examples/staticStructural/hingedTFrame.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/hingedTFrame"))
writedlm("test/newTestDataGenerators/staticStructural/hingedTFrame/u1_beam1.txt", test_data_matrix(u1_beam1))
writedlm("test/newTestDataGenerators/staticStructural/hingedTFrame/u1_beam2.txt", test_data_matrix(u1_beam2))
writedlm("test/newTestDataGenerators/staticStructural/hingedTFrame/u3_beam1.txt", test_data_matrix(u3_beam1))
writedlm("test/newTestDataGenerators/staticStructural/hingedTFrame/u3_beam2.txt", test_data_matrix(u3_beam2))
writedlm("test/newTestDataGenerators/staticStructural/hingedTFrame/F3_beam1.txt", test_data_matrix(F3_beam1))
writedlm("test/newTestDataGenerators/staticStructural/hingedTFrame/M2_beam1.txt", test_data_matrix(M2_beam1))

# Static analysis of the Lee frame (a right-angled frame) with a dead load
include("../examples/staticStructural/LeeFrameDeadLoad.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/LeeFrameDeadLoad"))
writedlm("test/newTestDataGenerators/staticStructural/LeeFrameDeadLoad/u1_atForce.txt", test_data_matrix(u1_atForce))
writedlm("test/newTestDataGenerators/staticStructural/LeeFrameDeadLoad/u3_atForce.txt", test_data_matrix(u3_atForce))

# Static analysis of the Lee frame (a right-angled frame) with a follower load
include("../examples/staticStructural/LeeFrameFollowerLoad.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/LeeFrameFollowerLoad"))
writedlm("test/newTestDataGenerators/staticStructural/LeeFrameFollowerLoad/u1_atForce.txt", test_data_matrix(u1_atForce))
writedlm("test/newTestDataGenerators/staticStructural/LeeFrameFollowerLoad/u3_atForce.txt", test_data_matrix(u3_atForce))

# Static analysis of two cantilevers connected at their tips by nonlinear doubly-attached springs in several directions, with and without rigid-body rotations of the assembly
include("../examples/staticStructural/nonlinearDoublyAttachedSprings.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/nonlinearDoublyAttachedSprings"))
writedlm("test/newTestDataGenerators/staticStructural/nonlinearDoublyAttachedSprings/tipDisp.txt", test_data_matrix(tipDisp))
writedlm("test/newTestDataGenerators/staticStructural/nonlinearDoublyAttachedSprings/clampLoads.txt", test_data_matrix(clampLoads))

# Static analysis of the pure bending test of the Pazy wing
include("../examples/staticStructural/PazyWingBendingTest.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/PazyWingBendingTest"))
writedlm("test/newTestDataGenerators/staticStructural/PazyWingBendingTest/tip_OOP.txt", test_data_matrix(tip_OOP))

# Static analysis of the coupled torsion-bending test of the Pazy wing
include("../examples/staticStructural/PazyWingTorsionTest.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/PazyWingTorsionTest"))
writedlm("test/newTestDataGenerators/staticStructural/PazyWingTorsionTest/tip_OOP.txt", test_data_matrix(tip_OOP))
writedlm("test/newTestDataGenerators/staticStructural/PazyWingTorsionTest/tip_twist.txt", test_data_matrix(tip_twist))

# Static analysis of a mid-loaded arch pinned at one end and clamped at the other
include("../examples/staticStructural/pinnedClampedArch.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/pinnedClampedArch"))
writedlm("test/newTestDataGenerators/staticStructural/pinnedClampedArch/u1.txt", test_data_matrix(u1_atForce))
writedlm("test/newTestDataGenerators/staticStructural/pinnedClampedArch/u3.txt", test_data_matrix(u3_atForce))

# Static analysis of a right-angled frame under a tip transverse follower force
include("../examples/staticStructural/rightAngledFrame.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/rightAngledFrame"))
writedlm("test/newTestDataGenerators/staticStructural/rightAngledFrame/tip_u1.txt", test_data_matrix(tip_u1))
writedlm("test/newTestDataGenerators/staticStructural/rightAngledFrame/tip_u3.txt", test_data_matrix(tip_u3))
writedlm("test/newTestDataGenerators/staticStructural/rightAngledFrame/tip_angle.txt", test_data_matrix(tip_angle))

# Static analysis of a right-angled frame under a 'buckling' load
include("../examples/staticStructural/rightAngledFrameBuckling.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/rightAngledFrameBuckling"))
writedlm("test/newTestDataGenerators/staticStructural/rightAngledFrameBuckling/tip_u2.txt", test_data_matrix(tip_u2))

# Static analysis of an L-frame with a doubly-attached spring and tip load
include("../examples/staticStructural/springedLFrame.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/springedLFrame"))
writedlm("test/newTestDataGenerators/staticStructural/springedLFrame/u3_b.txt", test_data_matrix(u3_b))

# Static analysis of a swept-back clamped beam with a driven flared folding wingtip
include("../examples/staticStructural/sweptBackFFWTWing.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/sweptBackFFWTWing"))
writedlm("test/newTestDataGenerators/staticStructural/sweptBackFFWTWing/phi.txt", test_data_matrix(ϕHinge))
writedlm("test/newTestDataGenerators/staticStructural/sweptBackFFWTWing/hingeMoment.txt", test_data_matrix(hingeMoment))

# Static analysis of a semi-circular arch with tangential follower force
include("../examples/staticStructural/tangentiallyForcedArch.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/tangentiallyForcedArch"))
writedlm("test/newTestDataGenerators/staticStructural/tangentiallyForcedArch/tip_u1.txt", test_data_matrix(tip_u1))
writedlm("test/newTestDataGenerators/staticStructural/tangentiallyForcedArch/tip_u3.txt", test_data_matrix(tip_u3))
writedlm("test/newTestDataGenerators/staticStructural/tangentiallyForcedArch/tip_angle.txt", test_data_matrix(tip_angle))

# Static analysis of a cantilever with tip follower transverse force 
include("../examples/staticStructural/tipFollowerForceCantilever.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/tipFollowerForceCantilever"))
writedlm("test/newTestDataGenerators/staticStructural/tipFollowerForceCantilever/tip_u1.txt", test_data_matrix(tip_u1))
writedlm("test/newTestDataGenerators/staticStructural/tipFollowerForceCantilever/tip_u3.txt", test_data_matrix(tip_u3))

# Static analysis of a cantilever with tip follower transverse force (force split over 2 BCs)
include("../examples/staticStructural/tipFollowerForceCantilever2.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/tipFollowerForceCantilever2"))
writedlm("test/newTestDataGenerators/staticStructural/tipFollowerForceCantilever2/tip_u1.txt", test_data_matrix(tip_u1))
writedlm("test/newTestDataGenerators/staticStructural/tipFollowerForceCantilever2/tip_u3.txt", test_data_matrix(tip_u3))

# Static analysis of a cantilever with tip moment
include("../examples/staticStructural/tipMomentCantilever.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/tipMomentCantilever"))
writedlm("test/newTestDataGenerators/staticStructural/tipMomentCantilever/tip_u1.txt", test_data_matrix(tip_u1))
writedlm("test/newTestDataGenerators/staticStructural/tipMomentCantilever/tip_u3.txt", test_data_matrix(tip_u3))
writedlm("test/newTestDataGenerators/staticStructural/tipMomentCantilever/tip_angle.txt", test_data_matrix(tip_angle))

# Static analysis of two cantilevers connected at their tips by a transverse doubly-attached spring, at several locations of the spring, against the analytical solution
include("../examples/staticStructural/transverseDoublyAttachedSpringCantilevers.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/transverseDoublyAttachedSpringCantilevers"))
writedlm("test/newTestDataGenerators/staticStructural/transverseDoublyAttachedSpringCantilevers/springDisp.txt", test_data_matrix(springDisp))
writedlm("test/newTestDataGenerators/staticStructural/transverseDoublyAttachedSpringCantilevers/tipDisp.txt", test_data_matrix(tipsDisp))
writedlm("test/newTestDataGenerators/staticStructural/transverseDoublyAttachedSpringCantilevers/tipTwist.txt", test_data_matrix(tipsTwist))

# Static analysis of a semi-circular arch with transverse follower force
include("../examples/staticStructural/transverselyForcedArch.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/transverselyForcedArch"))
writedlm("test/newTestDataGenerators/staticStructural/transverselyForcedArch/tip_u1.txt", test_data_matrix(tip_u1))
writedlm("test/newTestDataGenerators/staticStructural/transverselyForcedArch/tip_u3.txt", test_data_matrix(tip_u3))
writedlm("test/newTestDataGenerators/staticStructural/transverselyForcedArch/tip_angle.txt", test_data_matrix(tip_angle))

# Static analysis of a beam with triangular distributed load
include("../examples/staticStructural/triangleLoadBeam.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/triangleLoadBeam"))
writedlm("test/newTestDataGenerators/staticStructural/triangleLoadBeam/u3.txt", test_data_matrix(u3))
writedlm("test/newTestDataGenerators/staticStructural/triangleLoadBeam/F3.txt", test_data_matrix(F3))
writedlm("test/newTestDataGenerators/staticStructural/triangleLoadBeam/M2.txt", test_data_matrix(M2))

# Static analysis of a cantilever beam with a torsional spring attached between its middle and tip nodes, subjected to a torsional tip moment
include("../examples/staticStructural/twistDoublyAttachedSpringCantilever.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/twistDoublyAttachedSpringCantilever"))
writedlm("test/newTestDataGenerators/staticStructural/twistDoublyAttachedSpringCantilever/p1.txt", test_data_matrix(p1))
writedlm("test/newTestDataGenerators/staticStructural/twistDoublyAttachedSpringCantilever/M1.txt", test_data_matrix(M1))

# Static analysis of a beam with a central hinge driven by opposite end moments
include("../examples/staticStructural/drivenSpringedHingedBeam2.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/drivenSpringedHingedBeam2"))
writedlm("test/newTestDataGenerators/staticStructural/drivenSpringedHingedBeam2/u3.txt", test_data_matrix(u3))
writedlm("test/newTestDataGenerators/staticStructural/drivenSpringedHingedBeam2/p2.txt", test_data_matrix(p2))
writedlm("test/newTestDataGenerators/staticStructural/drivenSpringedHingedBeam2/M2.txt", test_data_matrix(M2))

# Static analysis of a swept cantilever under tip transverse loading
include("../examples/staticStructural/sweptBeamStaticLoading.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/staticStructural/sweptBeamStaticLoading"))
writedlm("test/newTestDataGenerators/staticStructural/sweptBeamStaticLoading/tip_u3.txt", test_data_matrix(tip_u3))
writedlm("test/newTestDataGenerators/staticStructural/sweptBeamStaticLoading/tip_twist.txt", test_data_matrix(tip_twist))