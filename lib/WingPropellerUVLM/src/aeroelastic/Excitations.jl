"""
    smooth_hann_pulse(time; start_time, duration)

Return a unit-amplitude Hann pulse that is zero outside the specified time
window and smoothly rises from and returns to zero inside it.
"""
function smooth_hann_pulse(time::Real; start_time::Real, duration::Real)
    start_time >= 0 || throw(ArgumentError("start_time must be nonnegative"))
    duration > 0 || throw(ArgumentError("duration must be positive"))
    start_time <= time <= start_time + duration || return 0.0
    phase = (time - start_time) / duration
    return 0.5 * (1.0 - cos(2.0 * pi * phase))
end

"""Create a generalized-load vector with the same Hann pulse on selected DOFs."""
function smooth_hann_pulse_load(time::Real, number_of_dofs::Integer,
    loaded_dofs::AbstractVector{<:Integer}; magnitude::Real,
    start_time::Real, duration::Real)

    number_of_dofs >= 0 || throw(ArgumentError("number_of_dofs must be nonnegative"))
    all(1 .<= loaded_dofs .<= number_of_dofs) ||
        throw(BoundsError(1:number_of_dofs, loaded_dofs))
    load = zeros(promote_type(Float64, typeof(magnitude)), number_of_dofs)
    load[loaded_dofs] .= magnitude * smooth_hann_pulse(
        time;
        start_time,
        duration,
    )
    return load
end
