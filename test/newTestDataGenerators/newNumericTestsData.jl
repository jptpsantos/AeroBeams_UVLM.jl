# Invariance under rigid-body motions and equilibrium of the nodal loads of a doubly-attached spring, for random nodal states and several locations of the spring
include("../examples/numeric/doublyAttachedSpringRigidBodyInvariance.jl")
# Residuals are checked against thresholds only (round-off level, random inputs), so no reference data is written

# Analytical Jacobian and inertia matrices of hinged beams, with and without aerodynamic surface, against finite differences
include("../examples/numeric/hingedBeamJacobianInertia.jl")
# Mismatch counts and inertia errors are checked against thresholds only (round-off level), so no reference data is written
