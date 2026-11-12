module ChangAeroBeams

using LinearAlgebra
using SparseArrays
using DelimitedFiles

# Resolve the two checked-out packages in the coupled-examples environment.
const REPOSITORY_ROOT = rstrip(normpath(joinpath(@__DIR__,"..","..","..","..","..")),('/','\\'))
const UVLM_ROOT = rstrip(normpath(joinpath(@__DIR__,"..","..","..")),('/','\\'))
const COUPLED_ENVIRONMENT = rstrip(normpath(joinpath(@__DIR__,"..","..")),('/','\\'))
COUPLED_ENVIRONMENT in LOAD_PATH || pushfirst!(LOAD_PATH,COUPLED_ENVIRONMENT)
import AeroBeams as AB
import WingPropellerUVLM as UVLM
import ForwardDiff
include(joinpath(@__DIR__,"..","..","AeroBeamsCoupling.jl"))
using .AeroBeamsCoupling

include(joinpath(@__DIR__,"..","..","chang_linear_aeroelastic","src","ChangAeroelastic.jl"))
const Chang = ChangAeroelastic

include("AeroBeamsWing.jl")
include("BeamPylon.jl")
include("ReducedRotor.jl")
include("TipRotor.jl")
include("AerodynamicInterface.jl")
include("CoupledSolver.jl")
include("Studies.jl")

export chang_aerobeams_defaults,build_chang_aerobeams,run_chang_aerobeams,
    solve_chang_aerobeams!,save_chang_aerobeams,initialize_mean_equilibrium!,
    ReducedRotor,TipRotor,rotor_mechanics,coupled_grids,aerodynamic_loads,
    build_wing,wing_node_poses,wing_node_parameters,wing_node_motion,
    run_chang_aerobeams_speed_sweep,perturbation_growth_rate

end
