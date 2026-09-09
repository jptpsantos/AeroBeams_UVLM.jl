# Mutable geometry, aerodynamic loads, and wake storage owned by one run.
# Constructing a workspace allocates fresh arrays; nothing is bound globally.

"""Allocate a fresh UVLM system and geometry/load buffers for `model`."""
function build_chang_workspace(model)
    uvlm = initialize_chang_uvlm(model.parameters, model.aerodynamic_options)
    count = model.parameters.Npropellers
    return merge(uvlm, (;
        T_hub_A_current = Vector{SVector{3,Float64}}(undef, count),
        T_load_A_current = Vector{SVector{3,Float64}}(undef, count),
        pitch_axis_A_current = Vector{SVector{3,Float64}}(undef, count),
        yaw_axis_A_current = Vector{SVector{3,Float64}}(undef, count),
    ))
end
