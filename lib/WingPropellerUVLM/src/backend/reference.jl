# --- reference quantities --- #

"""
    Reference(S, c, b, r, V[, rho])

Reference quantities.

**Arguments**
 - `S`: reference area
 - `c`: reference chord
 - `b`: reference span
 - `r`: reference location for all rotations/moments
 - `V`: reference velocity (magnitude)
 - `rho`: fluid density; defaults to `1.0` for backward compatibility
"""
struct Reference{TF}
    S::TF
    c::TF
    b::TF
    r::SVector{3, TF}
    V::TF
    rho::TF
end

function Reference(S, c, b, r, V, rho=1.0)
    TF = promote_type(typeof(S), typeof(c), typeof(b), eltype(r), typeof(V), typeof(rho))
    return Reference{TF}(S, c, b, r, V, rho)
end

Base.eltype(::Type{Reference{TF}}) where TF = TF
Base.eltype(::Reference{TF}) where TF = TF

Reference{TF}(r::Reference) where TF = Reference{TF}(r.S, r.c, r.b, r.r, r.V, r.rho)
Base.convert(::Type{Reference{TF}}, r::Reference) where {TF} = Reference{TF}(r)

# --- reference frames --- #

"""
    AbstractFrame

Supertype for the different possible reference frames used by this package.
"""
abstract type AbstractFrame end

"""
   Body <: AbstractFrame

Reference frame aligned with the global X-Y-Z axes
"""
struct Body <: AbstractFrame end

"""
    Stability <: AbstractFrame

Reference frame rotated from the body frame about the y-axis to be aligned with
the freestream `alpha`.
"""
struct Stability <: AbstractFrame end

"""
    Wind <: AbstractFrame

Reference frame rotated to be aligned with the freestream `alpha` and `beta`
"""
struct Wind <: AbstractFrame end

