"""Resolved SI flow and prescribed rotor operating point. Angles are radians."""
struct OperatingPoint
    airspeed::Float64
    density::Float64
    angleOfAttack::Float64
    sideslip::Float64
    angularSpeed::Float64
    referenceAirspeed::Float64
    referenceRotationRPM::Float64
    rotorSpeedPolicy::Symbol
end

function create_OperatingPoint(; airspeed, density=1.225, angleOfAttack=0.0,
    sideslip=0.0, rotationRPM=nothing, referenceAirspeed=airspeed,
    referenceRotationRPM=nothing, rotorSpeedPolicy=:prescribed)
    all(x -> x isa Real && !(x isa Bool) && isfinite(x),
        (airspeed,density,angleOfAttack,sideslip,referenceAirspeed)) ||
        throw(ArgumentError("Operating point must contain finite real values"))
    airspeed>0 && density>0 && referenceAirspeed>0 ||
        throw(ArgumentError("Airspeed, reference airspeed and density must be positive"))
    rotorSpeedPolicy in (:prescribed,:constantAdvanceRatio) ||
        throw(ArgumentError("rotorSpeedPolicy must be :prescribed or :constantAdvanceRatio"))
    if rotorSpeedPolicy==:constantAdvanceRatio
        isnothing(rotationRPM) || throw(ArgumentError("Specify referenceRotationRPM, not rotationRPM, for proportional scaling"))
        isnothing(referenceRotationRPM) && throw(ArgumentError("referenceRotationRPM is required"))
        rpm=referenceRotationRPM*airspeed/referenceAirspeed
    else
        isnothing(referenceRotationRPM) || throw(ArgumentError("Use rotationRPM for :prescribed speed"))
        rpm=something(rotationRPM,0.0)
        referenceRotationRPM=rpm
    end
    all(x -> x isa Real && !(x isa Bool) && isfinite(x),(rpm,referenceRotationRPM)) ||
        throw(ArgumentError("Rotor speeds must be finite"))
    return OperatingPoint(airspeed,density,angleOfAttack,sideslip,2pi*rpm/60,
        referenceAirspeed,referenceRotationRPM,rotorSpeedPolicy)
end
freestream(p::OperatingPoint) = Freestream(p.airspeed,p.angleOfAttack,p.sideslip,zeros(3))

