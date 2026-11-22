## AEROSOLVER COMPOSITE TYPES
# ------------------------------------------------------------------------------
"""

    AD (Automatic Differentiation) DerivationMethod composite type

"""
struct AD end
export AD


"""

    FD (Finite Differences) DerivationMethod composite type

# Fields
- method::FiniteDifferenceMethod
"""
@with_kw struct FD

    method::FiniteDifferences.FiniteDifferenceMethod = central_fdm(3,1)

    function FD(method::FiniteDifferences.FiniteDifferenceMethod)
        new(method)
    end

end
export FD
FD(FDFun::Function; nPoints::Int,adapt::Int=1) = begin

    @assert typeof(FDFun) in [typeof(FiniteDifferences.central_fdm),typeof(FiniteDifferences.forward_fdm),typeof(FiniteDifferences.backward_fdm)] "FDFun must be one of: central_fdm, forward_fdm, backward_fdm"
    @assert nPoints > 1 "nPoints must be greater than 1"
    @assert adapt >= 1 "adapt must be greater than or equal to 1"

    method = FDFun(nPoints,1,adapt=adapt)

    return FD(method)

end


"""

    QuasiSteady AeroSolver composite type

Quasi-steady aerodynamic solver for attached flow.

"""
struct QuasiSteady
    
    name::String
    nStates::Int
    nCirculatoryStates::Int
    nInertialStates::Int
    incompressibleInertialLoads::Bool

    function QuasiSteady()

        # Name of the aero solver
        name = "QuasiSteady"

        # Number of states per element
        nStates = nCirculatoryStates = nInertialStates = 0

        # Flag for incompressible inertial loads
        incompressibleInertialLoads = true

        return new(
            name,
            nStates,
            nCirculatoryStates,
            nInertialStates,
            incompressibleInertialLoads
        )

    end

end
export QuasiSteady


"""

    Indicial AeroSolver composite type

Unsteady aerodynamic solver for attached flow based on the indicial response of the airfoil.

# Keyword arguments
- `circulatoryIndicialFunction::Symbol`: indicial function of the circulatory loads, which can be `:Wagner` (default), `:Beddoes` or `:Jose`
- `inertialIndicialFunction::Symbol`: indicial function of the inertial (non-circulatory) loads (`:Leishman`)
- `incompressibleInertialLoads::Bool`: flag to compute the inertial loads with the incompressible (closed-form) expressions, which has no inertial states (default: `true`)
- `finiteWingLag::Bool`: flag to apply a finite-wing correction to the lag of the circulatory loads of the elements with tip correction (default: `true`). The tip loss factor ϖ models the steady loss of lift due to three-dimensional effects, but when it scales the whole circulatory response it also reduces the initial (impulsive) circulatory load of the Wagner-type function, 1-ΣAₙ = 1/2, which is that of the two-dimensional wing regardless of the aspect ratio (the three-dimensional effects take time to develop). With this option, the amplitudes Aₙ of the circulatory indicial function are scaled by κ = (1-(1-ΣAₙ)/ϖ)/ΣAₙ = 2 - 1/ϖ (not less than zero), so that the initial circulatory load, ϖ(1-κΣAₙ) = 1-ΣAₙ, is the two-dimensional one, while the steady one remains ϖ. The circulatory loads at the reduced frequencies of the flutter of wings of low aspect ratio are then closer to those of a three-dimensional unsteady potential flow solution, with less lag where the loss of lift is higher. Elements without tip correction (ϖ=1) are not affected
"""
struct Indicial{Nc,Ni}

    name::String
    circulatoryIndicialFunction::Symbol
    inertialIndicialFunction::Symbol
    incompressibleInertialLoads::Bool
    finiteWingLag::Bool
    nCirculatoryStates::Int
    nInertialStates::Int
    nStates::Int
    AC::SVector{Nc,Float64}
    bC::SVector{Nc,Float64}
    AI::SVector{Ni,Float64}
    bI::SVector{Ni,Float64}

    function Indicial(; 
        circulatoryIndicialFunction::Symbol=:Wagner,
        inertialIndicialFunction::Symbol=:Leishman, 
        incompressibleInertialLoads::Bool=true,
        finiteWingLag::Bool=true
        )

        # Validate
        @assert circulatoryIndicialFunction in [:Beddoes,:Wagner,:Jose] "circulatoryIndicialFunction must be one of: `:Beddoes`, `:Wagner`, `:Jose`"
        @assert inertialIndicialFunction in [:Leishman] "inertialIndicialFunction must be one of: `:Leishman`"

        # Name of the aero solver
        name = "Indicial"

        # Circulatory indicial parameters
        AC,bC = circulatory_indicial_parameters(circulatoryIndicialFunction)

        # Inertial indicial parameters
        AI,bI = inertial_indicial_parameters(inertialIndicialFunction)

        # Number of states per element
        nCirculatoryStates = length(AC)
        nInertialStates = incompressibleInertialLoads ?
                          0 : 
                          8
        nStates = nCirculatoryStates + nInertialStates

        # Set types
        Nc = length(AC)
        Ni = length(AI)

        return new{Nc,Ni}(
            name,
            circulatoryIndicialFunction,
            inertialIndicialFunction,
            incompressibleInertialLoads,
            finiteWingLag,
            nCirculatoryStates,
            nInertialStates,
            nStates,
            AC,
            bC,
            AI,
            bI
        )

    end

end
export Indicial


"""

    Inflow AeroSolver composite type

Unsteady aerodynamic solver based on the Peters's et al. inflow theory.

# Keyword arguments:
- `nInflowStates::Int`: number of inflow states (default: 6, can be either 4, 6, or 8)
- `circulatoryIndicialFunction::Symbol`: indicial function of the circulatory loads, which can be `:Wagner` (default), `:Beddoes` or `:Jose`
- `inertialIndicialFunction::Symbol`: indicial function of the inertial (non-circulatory) loads (`:Leishman`)
- `incompressibleInertialLoads::Bool`: flag to compute the inertial loads with the incompressible (closed-form) expressions, which has no inertial states (default: `true`)
- `finiteWingLag::Bool`: whether to include finite wing lag effects (default: `true`)
"""
struct Inflow{N,Nc,Ni,NN}

    name::String
    circulatoryIndicialFunction::Symbol
    inertialIndicialFunction::Symbol
    incompressibleInertialLoads::Bool
    finiteWingLag::Bool
    nInflowStates::Int
    nInertialStates::Int
    nStates::Int
    Aₚ::SMatrix{N,N,Float64,NN}
    AₚInv::SMatrix{N,N,Float64,NN}
    AₚInvVec::SVector{NN,Float64}
    AₚInvcₚ::SVector{N,Float64}
    bₚ::SVector{N,Float64}
    AC::SVector{Nc,Float64}
    bC::SVector{Nc,Float64}
    AI::SVector{Ni,Float64}
    bI::SVector{Ni,Float64}

    function Inflow(; 
        nInflowStates::Int=6,
        circulatoryIndicialFunction::Symbol=:Wagner,
        inertialIndicialFunction::Symbol=:Leishman,
        incompressibleInertialLoads::Bool=true,
        finiteWingLag::Bool=true,
        )

        # Validate
        @assert nInflowStates in [4,6,8] "select nInflowStates as 4, 6 or 8"
        @assert circulatoryIndicialFunction in [:Beddoes,:Wagner,:Jose] "circulatoryIndicialFunction must be one of: `:Beddoes`, `:Wagner`, `:Jose`"
        @assert inertialIndicialFunction in [:Leishman] "inertialIndicialFunction must be one of: `:Leishman`"

        # Name of the aero solver
        name = "Inflow"

        # Inflow arrays
        Aₚ,AₚInv,bₚ,cₚ = inflow_arrays(nInflowStates)
        AₚInvcₚ = AₚInv*cₚ
        AₚInvVec = reduce(vcat, AₚInv)

        # Number of states per element
        nInertialStates = incompressibleInertialLoads ?
                            0 :
                            8
        nStates = nInflowStates + nInertialStates

        # Circulatory indicial parameters (needed only to compute compressible inertial response time scales)
        AC,bC = circulatory_indicial_parameters(circulatoryIndicialFunction)

        # Inertial indicial parameters
        AI,bI = inertial_indicial_parameters(inertialIndicialFunction)

        # Set types
        N = nInflowStates
        Nc = length(AC)
        Ni = length(AI)

        return new{N,Nc,Ni,N*N}(
            name,
            circulatoryIndicialFunction,
            inertialIndicialFunction,
            incompressibleInertialLoads,
            finiteWingLag,
            nInflowStates,
            nInertialStates,
            nStates,
            Aₚ,
            AₚInv,
            AₚInvVec,
            AₚInvcₚ,
            bₚ,
            AC,
            bC,
            AI,
            bI
        )
    end

end
export Inflow


"""

    Incompressible modified Beddoes-Leishman AeroSolver composite type

Modified Beddoes-Leishman model matched to low-speed aerodynamics.

# Keyword arguments
- `circulatoryIndicialFunction::Symbol`: indicial function of the circulatory loads, which can be `:Wagner`, `:Beddoes` (default) or `:Jose`
- `inertialIndicialFunction::Symbol`: indicial function of the inertial (non-circulatory) loads (`:Leishman`)
- `incompressibleInertialLoads::Bool`: flag to compute the inertial loads with the incompressible (closed-form) expressions, which has no inertial states (default: `true`)
- `finiteWingLag::Bool`: whether to include finite wing lag effects (default: `true`)
- `qTransition::Real`: non-dimensional pitch rate for smoothing the transition between upstroke and downstroke. The transition between the upstroke and downstroke phases is smoothed over the non-dimensional pitch rate `qTransition`. A hard switch is recovered with `qTransition=0` (default: `1e-3`).
"""
struct BLi{Nc,Ni}

    name::String
    circulatoryIndicialFunction::Symbol
    inertialIndicialFunction::Symbol
    incompressibleInertialLoads::Bool
    finiteWingLag::Bool
    nCirculatoryStates::Int
    nInertialStates::Int
    nNonlinearStates::Int
    nStates::Int
    AC::SVector{Nc,Float64}
    bC::SVector{Nc,Float64}
    AI::SVector{Ni,Float64}
    bI::SVector{Ni,Float64}
    qTransition::Float64

    function BLi(;
        circulatoryIndicialFunction::Symbol=:Beddoes,
        inertialIndicialFunction::Symbol=:Leishman,
        incompressibleInertialLoads::Bool=true,
        finiteWingLag::Bool=true,
        qTransition::Real=1e-3
        )

        # Validate
        @assert circulatoryIndicialFunction in [:Beddoes,:Wagner,:Jose] "circulatoryIndicialFunction must be one of: `:Beddoes`, `:Wagner`, `:Jose`"
        @assert inertialIndicialFunction in [:Leishman] "inertialIndicialFunction must be one of: `:Leishman`"
        @assert qTransition >= 0 "qTransition must be non-negative"

        # Name of the aero solver
        name = "BLi"

        # Circulatory indicial parameters
        AC,bC = circulatory_indicial_parameters(circulatoryIndicialFunction)

        # Inertial indicial parameters
        AI,bI = inertial_indicial_parameters(inertialIndicialFunction)

        # Number of states per element
        nCirculatoryStates = length(AC)
        nInertialStates = incompressibleInertialLoads ? 
                          0 : 
                          8
        nNonlinearStates = 7
        nStates = nCirculatoryStates + nInertialStates + nNonlinearStates

        # Set types
        Nc = length(AC)
        Ni = length(AI)

        return new{Nc,Ni}(
            name,
            circulatoryIndicialFunction,
            inertialIndicialFunction,
            incompressibleInertialLoads,
            finiteWingLag,
            nCirculatoryStates,
            nInertialStates,
            nNonlinearStates,
            nStates,
            AC,
            bC,
            AI,
            bI,
            qTransition
        )

    end

end
export BLi


"""

    Original Beddoes-Leishman AeroSolver composite type

Interpretation of the original Beddoes-Leishman model for dynamic stall.

# Keyword arguments
- `circulatoryIndicialFunction::Symbol`: indicial function of the circulatory loads, which can be `:Wagner`, `:Beddoes` (default) or `:Jose`
- `inertialIndicialFunction::Symbol`: indicial function of the inertial (non-circulatory) loads (`:Leishman`)
- `incompressibleInertialLoads::Bool`: flag to compute the inertial loads with the incompressible (closed-form) expressions, which has no inertial states (default: `true`)
- `finiteWingLag::Bool`: whether to include finite wing lag effects (default: `true`)
- `qTransition::Real`: non-dimensional pitch rate for smoothing the transition between upstroke and downstroke. The transition between the upstroke and downstroke phases is smoothed over the non-dimensional pitch rate `qTransition`. A hard switch is recovered with `qTransition=0` (default: `1e-3`).
"""
struct BLo{Nc,Ni}

    name::String
    circulatoryIndicialFunction::Symbol
    inertialIndicialFunction::Symbol
    incompressibleInertialLoads::Bool
    finiteWingLag::Bool
    nCirculatoryStates::Int
    nInertialStates::Int
    nNonlinearStates::Int
    nStates::Int
    AC::SVector{Nc,Float64}
    bC::SVector{Nc,Float64}
    AI::SVector{Ni,Float64}
    bI::SVector{Ni,Float64}
    qTransition::Float64

    function BLo(;
        circulatoryIndicialFunction::Symbol=:Beddoes,
        inertialIndicialFunction::Symbol=:Leishman,
        incompressibleInertialLoads::Bool=false,
        finiteWingLag::Bool=true,
        qTransition::Real=1e-3
        )

        # Validate
        @assert circulatoryIndicialFunction in [:Beddoes,:Wagner,:Jose] "circulatoryIndicialFunction must be one of: `:Beddoes`, `:Wagner`, `:Jose`"
        @assert inertialIndicialFunction in [:Leishman] "inertialIndicialFunction must be one of: `:Leishman`"
        @assert qTransition >= 0 "qTransition must be non-negative"

        # Name of the aero solver
        name = "BLo"

        # Circulatory indicial parameters
        AC,bC = circulatory_indicial_parameters(circulatoryIndicialFunction)

        # Inertial indicial parameters
        AI,bI = inertial_indicial_parameters(inertialIndicialFunction)

        # Number of states per element
        nCirculatoryStates = length(AC)
        nInertialStates = incompressibleInertialLoads ? 
                          0 : 
                          8
        nNonlinearStates = 5
        nStates = nCirculatoryStates + nInertialStates + nNonlinearStates

        # Set types
        Nc = length(AC)
        Ni = length(AI)

        return new{Nc,Ni}(
            name,
            circulatoryIndicialFunction,
            inertialIndicialFunction,
            incompressibleInertialLoads,
            finiteWingLag,
            nCirculatoryStates,
            nInertialStates,
            nNonlinearStates,
            nStates,
            AC,
            bC,
            AI,
            bI,
            qTransition
        )

    end

end
export BLo


"""

    TableLookup FlapAeroSolver composite type

Flap-induced aerodynamic loads solver using table lookup.

"""
struct TableLookup

    nStates::Int

    function TableLookup()

        # Number of states per element
        nStates = 0

        return new(nStates)

    end

end
export TableLookup


"""

    ThinAirfoilTheory FlapAeroSolver composite type

Flap-induced aerodynamic loads solver using thin airfoil theory.

# Keyword arguments
- `circulatoryIndicialFunction::Symbol`: indicial function of the circulatory loads, which can be `:Wagner` (default), `:Beddoes` or `:Jose`
"""
struct ThinAirfoilTheory{N}

    circulatoryIndicialFunction::Symbol
    nStates::Int
    ACf::SVector{N,Float64}
    bCf::SVector{N,Float64}

    function ThinAirfoilTheory(circulatoryIndicialFunction::Symbol=:Wagner)

        # Circulatory indicial parameters
        ACf,bCf = circulatory_indicial_parameters(circulatoryIndicialFunction)

        # Number of states per element
        N = nStates = length(bCf)

        return new{N}(
            circulatoryIndicialFunction,
            nStates,
            ACf,
            bCf
        )

    end

end
export ThinAirfoilTheory


"""

    IndicialGust GustAeroSolver composite type

Gust-induced aerodynamic loads solver based on indicial theory.

# Keyword arguments
- `indicialFunction::Symbol`: indicial function of the gust circulatory loads, which can be `:QuasiSteady`, `:Kussner` (default) or `:BerciRighi`.
"""
struct IndicialGust{N}

    indicialFunction::Symbol
    nStates::Int
    AG::SVector{N,Float64}
    bG::SVector{N,Float64}
    AGbG::SVector{N,Float64}

    function IndicialGust(indicialFunction::Symbol=:Kussner)

        # Validate
        @assert indicialFunction in [:QuasiSteady,:Kussner,:BerciRighi] "indicialFunction must be one of: `:QuasiSteady`, `:Kussner`, `:BerciRighi`"

        # Set number of states and indicial parameters arrays according to indicial function
        AG,bG = gust_circulatory_indicial_parameters(indicialFunction)

        # Number of gust states
        N = nStates = length(bG)

        # Useful functions of the indicial parameters arrays
        AGbG = AG .* bG

        return new{N}(
            indicialFunction,
            nStates,
            AG,
            bG,
            AGbG
        )

    end

end
export IndicialGust


# Gets the circulatory indicial parameters given the indicial function
function circulatory_indicial_parameters(circulatoryIndicialFunction)
    ACbC = Dict(
        :Beddoes => (SVector(0.3, 0.7), SVector(0.14, 0.53)),
        :Wagner  => (SVector(0.165, 0.335), SVector(0.0455, 0.3)),
        :Jose    => (SVector(0.3493, 0.6507), SVector(0.0984, 0.7759))
    )
    return get(ACbC, circulatoryIndicialFunction, nothing)
end


# Gets the inertial indicial parameters given the indicial function
function inertial_indicial_parameters(inertialIndicialFunction)
    AIbI = Dict(
        :Leishman => (SVector(1.5, -0.5, 1.0), SVector(0.25, 0.1, 5.0)),
    )
    return get(AIbI, inertialIndicialFunction, nothing)
end


# Gets the gust circulatory indicial parameters given the indicial function
function gust_circulatory_indicial_parameters(indicialFunction)
    AGbG = Dict(
        :QuasiSteady   => (zeros(SVector{0,Float64}), zeros(SVector{0,Float64})),
        :Kussner       => (SVector(0.5, 0.5), SVector(0.13, 1.0)),
        :BerciRighi    => (SVector(0.3694, 0.0550, 0.2654, 0.1829, 0.0861, 0.0412), SVector(0.3733, 0.0179, 0.1096, 1.6003, 10.428, 170.93))
    )
    return get(AGbG, indicialFunction, nothing)
end


# Computes the fixed state arrays from Peters' inflow theory
function inflow_arrays(N::Int)

    # D matrix
    D = zeros(N,N)
    D[1,2] = -1/2
    D[end,end-1] = 1/(2*N)
    for n=2:N-1
        D[n,n+1] = -1/(2*n)
        D[n,n-1] = 1/(2*n)
    end

    # b vector
    b = zeros(N)
    for n=1:N-1
        F = factorial(N+n-1)/factorial(N-n-1)/(factorial(n))^2
        b[n] = (-1)^(n-1)*F
    end
    b[N] = (-1)^(N+1) 

    # c vector
    c = zeros(N)
    for n=1:N
        c[n] = 2/n
    end

    # d vector
    d = zeros(N)
    d[1] = 1/2

    # A matrix and its inverse
    A = D + d*b' + c*d' + 1/2*c*b'
    AInv = inv(A)

    # Set as static arrays
    A = SMatrix{N,N,Float64,N*N}(A)
    AInv = SMatrix{N,N,Float64,N*N}(AInv)
    b = SVector{N,Float64}(b)
    c = SVector{N,Float64}(c)

    return A,AInv,b,c

end

# Empty solvers
struct NoAeroSolver end
struct NoFlapAeroSolver end
struct NoGustAeroSolver end

# Type aliases for solver subtypes
const AeroSolver = Union{QuasiSteady, Indicial, Inflow, BLi, BLo, NoAeroSolver}
const AttachedSolver = Union{QuasiSteady, Indicial, Inflow}
const UnsteadyAttachedSolver = Union{Indicial, Inflow}
const IndicialSolver = Union{Indicial, BLi, BLo}
const UnsteadySolver = Union{Indicial, Inflow, BLi, BLo}
const DynamicStallSolver = Union{BLi, BLo}

const FlapAeroSolver = Union{TableLookup, ThinAirfoilTheory, NoFlapAeroSolver}

const GustAeroSolver = Union{IndicialGust, NoGustAeroSolver}

const DerivationMethod = Union{AD,FD}

export AeroSolver, AttachedSolver, UnsteadyAttachedSolver, IndicialSolver, UnsteadySolver, DynamicStallSolver, FlapAeroSolver, GustAeroSolver, DerivationMethod

# Default gust solver
_defaultGustAeroSolver(solver) = solver isa QuasiSteady ? 
                                  IndicialGust(:QuasiSteady) : 
                                  IndicialGust(:Kussner)