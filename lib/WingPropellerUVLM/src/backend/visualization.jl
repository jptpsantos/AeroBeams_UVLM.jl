# TODO: Add grid based visualization
"""
    write_vtk(name, system; write_surfaces = true, write_wakes = false, kwargs...)

Write geometry from surfaces and/or wakes to Paraview files for visualization.

# Arguments
 - `name`: Base name for the generated files
 - `system`: System object containing surfaces and/or wakes

# Keyword Arguments:
    - `write_surfaces = true`: Flag indicating whether to write surface geometry
    - `write_wakes = false`: Flag indicating whether to write wake geometry
    - `symmetric`: (required if `surface_properties` is provided) Flags indicating whether a
        mirror image (across the X-Z plane) was used when calculating induced velocities
        for each surface.
    - `trailing_vortices`: Flag indicating whether the model uses trailing vortices.
        Defaults to `true` when wake panels are absent, `false` otherwise
    - `xhat`: Direction in which trailing vortices extend if used. Defaults to [1, 0, 0].
    - `wake_length`: Distance to extend trailing vortices. Defaults to 10
    - `metadata`: Dictionary of metadata to include in generated files
"""


function write_vtk(name::String,system::System; write_surfaces=true,write_wakes=false,kwargs...)

    if write_surfaces && write_wakes
        write_vtk(name,system.surfaces,system.wakes,system.properties; symmetric=system.symmetric,kwargs...)
    elseif write_surfaces
        write_vtk(name,system.surfaces,system.properties; symmetric=system.symmetric,kwargs...)
    elseif write_wakes
        write_vtk(name,system.wakes; symmetric=system.symmetric,kwargs...)
    else
        error("At least one of `write_surfaces` or `write_wakes` must be true")
    end

    return nothing
end

"""
    write_vtk(name, surfaces, [surface_properties]; kwargs...)
    write_vtk(name, wakes; kwargs...)
    write_vtk(name, surfaces, wakes, [surface_properties]; kwargs...)

Write geometry from surfaces and/or wakes to Paraview files for visualization.

# Arguments
 - `name`: Base name for the generated files
 - `surfaces`:
   - Vector of grids of shape (3, nc+1, ns+1) which represent lifting surfaces
   or
   - Vector of matrices of shape (nc, ns) containing surface panels (see
    [`SurfacePanel`](@ref))
   where `nc` is the number of chordwise panels and `ns` is the number of
   spanwise panels
 - `wakes`: (optional) Vector of wakes corresponding to each surface, represented
    by matrices of wake panels (see [`WakePanel`](@ref)) of shape (nw, ns) where
    `nw` is the number of chordwise wake panels and `ns` is the number of
    spanwise panels.
 - `surface_properties`: (optional) Vector of surface panel properties for each
    surface, stored as matrices of panel properties (see [`PanelProperties`](@ref))
    of shape (nc, ns) where `nc` is the number of chordwise panels and `ns` is
    the number of spanwise panels

# Keyword Arguments:
 - `symmetric`: (required if `surface_properties` is provided) Flags indicating whether a
    mirror image (across the X-Z plane) was used when calculating induced velocities
    for each surface.
 - `trailing_vortices`: Flag indicating whether the model uses trailing vortices.
    Defaults to `true` when wake panels are absent, `false` otherwise
 - `xhat`: Direction in which trailing vortices extend if used. Defaults to [1, 0, 0].
 - `wake_length`: Distance to extend trailing vortices. Defaults to 10
 - `metadata`: Dictionary of metadata to include in generated files
"""
write_vtk(name,surfaces::AbstractVector{<:AbstractMatrix},args...; kwargs...)

function write_vtk(name,surfaces::AbstractVector{<:AbstractMatrix{<:SurfacePanel}},
    properties=fill(nothing,length(surfaces)); symmetric=fill(nothing,length(surfaces)),kwargs...)

    # create paraview multiblock file
    vtk_multiblock(name) do vtmfile
        # loop through all surfaces
        for i = 1:length(surfaces)
            # add paraview files corresponding to the surface to the multiblock file
            write_vtk!(vtmfile,surfaces[i],properties[i]; symmetric=symmetric[i],kwargs...)
        end
    end

    return nothing
end

function write_vtk(name,wakes::AbstractVector{<:AbstractMatrix{<:WakePanel}}; kwargs...)

    # create paraview multiblock file
    vtk_multiblock(name) do vtmfile
        # loop through all wakes
        for i = 1:length(wakes)
            # add paraview files corresponding to the wake to the multiblock file
            write_vtk!(vtmfile,wakes[i]; kwargs...,surface_index=i)
        end
    end

    return nothing
end

function write_vtk(name,surfaces::AbstractVector{<:AbstractMatrix{<:SurfacePanel}},
    wakes::AbstractVector{<:AbstractMatrix{<:WakePanel}},properties=nothing; kwargs...)

    # create multiblock file
    vtk_multiblock(name) do vtmfile

        # loop through all surfaces
        for i = 1:length(surfaces)

            # extract circulation at the trailing edge of `surface`
            if isnothing(properties)
                surface_circulation = zeros(size(surfaces[i],2))
            else
                surface_circulation = getproperty.(properties[i][end,:],:gamma)
            end

            # extract circulation at the leading edge of `wake`
            wake_circulation = getproperty.(wakes[i][1,:],:gamma)

            # add paraview files corresponding to the surface to the multiblock file
            write_vtk!(vtmfile,surfaces[i],properties[i]; wake_circulation,
                trailing_edge=isempty(wakes[i]),kwargs...,trailing_vortices=false)

            # add paraview files corresponding to the wake to the multiblock file
            write_vtk!(vtmfile,wakes[i]; surface_circulation,kwargs...)
        end
    end

    return nothing
end

"""
    write_vtk(name, surface_history, property_history, wake_history; kwargs...)

Writes unsteady simulation geometry to Paraview files for visualization.

# Arguments
 - `name`: Base name for the generated files
 - `surface_history`: Vector of surfaces at each time step, where each surface is
    represented by a matrix of surface panels (see [`SurfacePanel`](@ref)) of shape
    (nc, ns) where `nc` is the number of chordwise panels and `ns` is the number
    of spanwise panels
 - `property_history`: Vector of surface properties for each surface at each
    time step, where surface properties are represented by a matrix of panel
    properties (see [`PanelProperties`](@ref)) of shape (nc, ns) where `nc` is
    the number of chordwise panels and `ns` is the number of spanwise panels
 - `wake_history`: Vector of wakes corresponding to each surface at each time step,
    where each wake is represented by a matrix of wake panels (see [`WakePanel`](@ref))
    of shape (nw, ns) where `nw` is the number of chordwise wake panels and
    `ns` is the number of spanwise panels.
 - `dt`: Time step vector

# Keyword Arguments:
 - `symmetric`: (required if `properties` is provided) Flags indicating whether a
    mirror image (across the X-Z plane) was used when calculating induced velocities
    for each surface.
 - `wake_length`: Distance to extend trailing vortices. Defaults to 10
 - `metadata`: Dictionary of metadata to include in generated files
"""
function write_vtk(name,surface_history::AbstractVector{<:AbstractVector{<:AbstractMatrix}},
    property_history::AbstractVector{<:AbstractVector{<:AbstractMatrix}},
    wake_history::AbstractVector{<:AbstractVector{<:AbstractMatrix}},dt;
    symmetric=fill(nothing,length(surface_history[1])),kwargs...)

    symmetric = isa(symmetric,Number) ? fill(symmetric,length(surface_history[1])) : symmetric

    # create paraview collection file
    paraview_collection(name) do pvdfile

        # construct time vector
        time = cumsum(dt)

        # loop through each time step
        for it = 1:length(time)

            # construct multiblock file for each time step
            vtk_multiblock(name * "-step$it") do vtmfile

                # loop through all surfaces
                for i = 1:length(surface_history[it])

                    # extract circulation at the trailing edge of `surface`
                    surface_circulation = getproperty.(property_history[it][i][end,:],:gamma)

                    # extract circulation at the leading edge of `wake`
                    if isempty(wake_history[it][i])
                        wake_circulation = zeros(size(surface_history[it][i],2))
                    else
                        wake_circulation = getproperty.(wake_history[it][i][1,:],:gamma)
                    end

                    # add paraview files corresponding to the surface to the multiblock file
                    write_vtk!(vtmfile,surface_history[it][i],property_history[it][i];
                        wake_circulation,
                        trailing_edge=isempty(wake_history[it][i]),
                        symmetric=symmetric[i],
                        kwargs...,
                        trailing_vortices=false)

                    # add paraview files corresponding to the wake to the multiblock file
                    write_vtk!(vtmfile,wake_history[it][i];
                        surface_circulation,
                        symmetric=symmetric[i],
                        kwargs...)
                end

                # add multiblock file to the paraview collection file
                pvdfile[time[it]] = vtmfile
            end
        end
    end

    return nothing
end

"""
    write_vtk!(vtmfile, surface, [surface_properties]; kwargs...)

Writes geometry to Paraview files for visualization.

# Arguments
 - `vtmfile`: Multiblock file handle
 - `surface`: Matrix of surface panels (see [`SurfacePanel`](@ref)) of shape
    (nc, ns) where `nc` is the number of chordwise panels and `ns` is the number
    of spanwise panels
 - `surface_properties`: (optional) Matrix of panel properties for each non-wake panel
    where each element of the matrix is of type [`PanelProperties`](@ref).

# Keyword Arguments:
 - `symmetric`: (required if `properties` is provided) Flag indicating whether a
    mirror image (across the X-Z plane) was used when calculating induced velocities.
 - `trailing_vortices = true`: Flag indicating whether the model uses trailing vortices
 - `xhat = [1, 0, 0]`: Direction in which trailing vortices extend if used
 - `wake_length = 10`: Distance to extend trailing vortices
 - `wake_circulation = zeros(size(surfaces, 2))`: Contribution to the trailing
    edge circulation from the wake attached to this surface
 - `metadata = Dict()`: Dictionary of metadata to include in generated files
"""
function write_vtk!(vtmfile,surface::AbstractMatrix{<:SurfacePanel},properties=nothing;
    symmetric=nothing,
    trailing_vortices=true,
    trailing_edge=true,
    xhat=SVector(1,0,0),
    wake_length=10,
    wake_circulation=zeros(size(surface,2)),
    metadata=Dict())

    # get float type
    TF = eltype(eltype(surface))

    # check to make sure `symmetric` is provided if `properties` is provided
    @assert !(!isnothing(properties) && isnothing(symmetric)) "Keyword argument `symmetric` is required when optional argument `properties` is provided"

    # grid dimensions
    nc,ns = size(surface)
    N = length(surface)

    # extract geometry as a grid
    xyz = Array{TF,4}(undef,3,nc + 1,ns + 1,1)
    for (i,I) in enumerate(CartesianIndices((nc + 1,ns + 1)))
        if I[1] <= nc && I[2] <= ns
            xyz[:,I,1] = top_left(surface[I[1],I[2]])
        elseif I[1] <= nc && I[2] == ns + 1
            xyz[:,I,1] = top_right(surface[I[1],I[2]-1])
        elseif I[1] == nc + 1 && I[2] <= ns
            xyz[:,I,1] = bottom_left(surface[I[1]-1,I[2]])
        else # I[1] == nc + 1 && I[2] == ns + 1
            xyz[:,I,1] = bottom_right(surface[I[1]-1,I[2]-1])
        end
    end

    # convert grid to points
    points = reshape(xyz,3,:)

    # now extract bound vortex geometries
    li = LinearIndices((size(xyz,2),size(xyz,3)))
    # horizontal bound vortices
    linesH = [MeshCell(PolyData.Lines(),[
        li[i,j],
        li[i,j+1]
    ]) for j = 1:ns for i = 1:nc]
    # vertical bound vortices
    linesV = [MeshCell(PolyData.Lines(),[
        li[i,j],
        li[i+1,j]
    ]) for j = 1:ns+1 for i = 1:nc]

    # now extract data (if applicable)
    if !isnothing(properties)
        # horizontal bound vortex circulation strength
        gammaH = Matrix{TF}(undef,nc,ns)
        for i = 1:nc,j = 1:ns
            previousGamma = i == 1 ? 0.0 : properties[i-1,j].gamma
            currentGamma = properties[i,j].gamma
            # check if we need to account for symmetry
            if symmetric && on_symmetry_plane(top_left(surface[i,j]),top_right(surface[i,j]))
                gammaH[i,j] = 0.0
            else
                gammaH[i,j] = currentGamma - previousGamma
            end
        end

        # horizontal bound vortex force coefficient
        cf_h = Array{TF}(undef,3,nc,ns)
        for i = 1:nc,j = 1:ns

            # check if we need to account for symmetry
            if symmetric && on_symmetry_plane(top_left(surface[i,j]),top_right(surface[i,j]))
                cf_h[:,i,j] .= 0.0
            else
                cf_h[:,i,j] = properties[i,j].cfb
            end
        end

        # horizontal bound vortex local velocity
        v_h = Array{TF}(undef,3,nc,ns)
        for i = 1:nc,j = 1:ns
            v_h[:,i,j] = properties[i,j].velocity
        end

        # vertical bound vortex circulation strength
        gammaV = Matrix{TF}(undef,nc,ns + 1)
        for i = 1:nc
            currentGamma = properties[i,1].gamma
            # check if we need to account for symmetry
            if symmetric && on_symmetry_plane(bottom_left(surface[i,1]),top_left(surface[i,1]))
                gammaV[i,1] = 0.0
            else
                gammaV[i,1] = -currentGamma
            end
            for j = 2:ns
                previousGamma = currentGamma
                currentGamma = properties[i,j].gamma
                # check if we need to account for symmetry
                if symmetric && on_symmetry_plane(bottom_left(surface[i,j]),top_left(surface[i,j]))
                    gammaV[i,j] = 0.0
                else
                    gammaV[i,j] = previousGamma - currentGamma
                end
            end
            previousGamma = currentGamma
            # check if we need to account for symmetry
            if symmetric && on_symmetry_plane(bottom_right(surface[i,end]),top_right(surface[i,end]))
                gammaV[i,end] = 0.0
            else
                gammaV[i,end] = previousGamma
            end
        end

        # vertical bound vortex force coefficient
        cf_v = Array{TF}(undef,3,nc,ns + 1)
        for i = 1:nc
            #  check if we need to account for symmetry
            if symmetric && on_symmetry_plane(bottom_left(surface[i,1]),top_left(surface[i,1]))
                cf_v[:,i,1] .= 0.0
            else
                cf_v[:,i,1] = properties[i,1].cfl
            end
            for j = 2:ns
                previousCf = properties[i,j-1].cfr
                currentCf = properties[i,j].cfl
                #  check if we need to account for symmetry
                if symmetric && on_symmetry_plane(bottom_left(surface[i,j]),top_left(surface[i,j]))
                    cf_v[i,j] = 0.0
                else
                    cf_v[:,i,j] = previousCf + currentCf
                end
            end
            #  check if we need to account for symmetry
            if symmetric && on_symmetry_plane(bottom_right(surface[i,end]),top_right(surface[i,end]))
                cf_v[:,i,end] .= 0.0
            else
                cf_v[:,i,end] = properties[i,end].cfr
            end
        end
    end

    # horizontal bound vortices
    vtk_grid(vtmfile,points,linesH) do vtkfile

        # add metadata
        for (key,value) in pairs(metadata)
            vtkfile[string(key)] = value
        end


        if !isnothing(properties)
            # circulation strength
            vtkfile["circulation"] = reshape(gammaH,:)

            # local velocity
            vtkfile["velocity"] = reshape(v_h,3,:)

            # force coefficient
            vtkfile["force"] = reshape(cf_h,3,:)
        end
    end

    # vertical bound vortices
    vtk_grid(vtmfile,points,linesV) do vtkfile

        # add metadata
        for (key,value) in pairs(metadata)
            vtkfile[string(key)] = value
        end

        if !isnothing(properties)
            # circulation strength
            vtkfile["circulation"] = reshape(gammaV,:)

            # force coefficient
            vtkfile["force"] = reshape(cf_v,3,:)
        end
    end

    # --- control points ---

    pointsCp = Matrix{TF}(undef,3,N)
    for i = 1:N
        ipoint = i
        pointsCp[:,ipoint] = controlpoint(surface[i])
    end

    cellsCp = [MeshCell(PolyData.Verts(),[i]) for i = 1:N]

    vtk_grid(vtmfile,pointsCp,cellsCp) do vtkfile

        # add metadata
        for (key,value) in pairs(metadata)
            vtkfile[string(key)] = value
        end

        # add normal
        data = Matrix{TF}(undef,3,N)
        for i = 1:length(surface)
            data[:,i] = normal(surface[i])
        end
        vtkfile["normal"] = data

    end

    # trailing edge and/or trailing vortices
    if trailing_vortices

        # trailing vortices points as a grid
        xyzT = Array{TF}(undef,3,2,ns + 1)
        for j = 1:ns
            xyzT[:,1,j] = bottom_left(surface[end,j])
            xyzT[:,2,j] = xyzT[:,1,j] + wake_length * xhat
        end
        xyzT[:,1,end] = bottom_right(surface[end,end])
        xyzT[:,2,end] = xyzT[:,1,end] + wake_length * xhat

        # generate points
        pointsT = reshape(xyzT,3,:)

        li = LinearIndices((2,ns + 1))

        # horizontal bound vortices
        linesH = [MeshCell(PolyData.Lines(),[li[1,j],li[1,j+1]]) for j = 1:ns]

        # vertical bound vortices
        linesV = [MeshCell(PolyData.Lines(),[li[1,j],li[2,j]]) for j = 1:ns+1]

        # combine
        linesT = vcat(linesH,linesV)

        # now extract data (if applicable)
        if !isnothing(properties)
            # horizontal bound vortex circulation strength
            gammaH = zeros(ns)

            # vertical bound vortex circulation strength
            gammaV = Vector{TF}(undef,ns + 1)
            currentGamma = properties[end,1].gamma

            # check if we need to account for symmetry
            if symmetric && on_symmetry_plane(bottom_left(surface[end,1]))
                gammaV[1] = 0.0
            else
                gammaV[1] = -currentGamma
            end

            for j = 2:ns
                previousGamma = currentGamma
                currentGamma = properties[end,j].gamma

                # check if we need to account for symmetry
                if symmetric && on_symmetry_plane(bottom_left(surface[end,j]))
                    gammaV[j] = 0.0
                else
                    gammaV[j] = previousGamma - currentGamma
                end

            end
            previousGamma = currentGamma

            # check if we need to account for symmetry
            if symmetric && on_symmetry_plane(bottom_right(surface[end,end]))
                gammaV[end] = 0.0
            else
                gammaV[end] = previousGamma
            end

            # combine
            gammaT = vcat(gammaH,gammaV)
        end

    else
        # only generate if we plan to include the trailing edge
        if trailing_edge

            # generate points
            pointsT = Matrix{TF}(undef,3,ns + 1)
            for j = 1:ns
                pointsT[:,j] = bottom_left(surface[end,j])
            end
            pointsT[:,end] = bottom_right(surface[end,end])

            # horizontal bound vortices
            linesT = [MeshCell(PolyData.Lines(),j:j+1) for j = 1:ns]

            # now extract data (if applicable)
            if !isnothing(properties)
                # horizontal bound vortex circulation strength
                gammaT = Vector{TF}(undef,ns)
                for j = 1:ns
                    # check if we need to account for symmetry
                    if symmetric && on_symmetry_plane(bottom_left(surface[end,j]),bottom_right(surface[end,j]))
                        gammaT[j] = 0.0
                    else
                        gammaT[j] = wake_circulation[j] - properties[end,j].gamma
                    end
                end
            end
        end

    end

    if trailing_vortices || trailing_edge

        # trailing edge and trailing vortices
        vtk_grid(vtmfile,pointsT,linesT) do vtkfile

            # add metadata
            for (key,value) in pairs(metadata)
                vtkfile[string(key)] = value
            end

            if !isnothing(properties)
                # circulation strength
                vtkfile["circulation"] = gammaT
            end
        end
    end

    return nothing
end

"""
    write_vtk!(vtmfile, wake; kwargs...)

Writes geometry to Paraview files for visualization.

# Arguments
 - `vtmfile`: Paraview file handle
 - `wake`: Matrix of wake panels (see [`WakePanel`](@ref)) of shape (nw, ns)
    where `nw` is the number of chordwise wake panels and `ns` is the number of
    spanwise panels

# Keyword Arguments:
 - `symmetric`: (required) Flag indicating whether a mirror image (across the
    X-Z plane) was used when calculating induced velocities.
 - `trailing_vortices = false`: Flag indicating whether the model uses trailing vortices
 - `xhat = [1, 0, 0]`: Direction in which trailing vortices extend if used
 - `wake_length = 10`: Distance to extend trailing vortices
 - `surface_circulation = zeros(size(wake, 2))`: Contribution to the leading edge
    circulation from the surface attached to this wake.
 - `metadata = Dict()`: Dictionary of metadata to include in generated files
"""
function write_vtk!(vtmfile,wake::AbstractMatrix{<:WakePanel};
    symmetric,
    trailing_vortices=false,
    xhat=SVector(1,0,0),
    wake_length=10,
    surface_circulation=zeros(size(wake,2)),
    metadata=Dict())

    # do nothing if no wake panels are present
    if isempty(wake)
        return vtmfile
    end

    # extract float type
    TF = eltype(eltype(wake))

    # get wake dimensions
    nc,ns = size(wake)
    N = length(wake)

    # extract geometry as a grid
    xyz = Array{TF,4}(undef,3,nc + 1,ns + 1,1)
    for (i,I) in enumerate(CartesianIndices((nc + 1,ns + 1)))
        if I[1] <= nc && I[2] <= ns
            xyz[:,I,1] = top_left(wake[I[1],I[2]])
        elseif I[1] <= nc && I[2] == ns + 1
            xyz[:,I,1] = top_right(wake[I[1],I[2]-1])
        elseif I[1] == nc + 1 && I[2] <= ns
            xyz[:,I,1] = bottom_left(wake[I[1]-1,I[2]])
        else # I[1] == nc + 1 && I[2] == ns + 1
            xyz[:,I,1] = bottom_right(wake[I[1]-1,I[2]-1])
        end
    end

    # convert to points
    points = reshape(xyz,3,:)

    # now extract bound vortex geometries
    li = LinearIndices(size(xyz)[2:3])
    # horizontal bound vortices
    linesH = [MeshCell(PolyData.Lines(),[
        li[i,j],
        li[i,j+1]
    ]) for j = 1:ns for i = 1:nc]
    # vertical bound vortices
    linesV = [MeshCell(PolyData.Lines(),[
        li[i,j],
        li[i+1,j]
    ]) for j = 1:ns+1 for i = 1:nc]

    # horizontal bound vortex circulation strength
    gammaH = Matrix{TF}(undef,nc,ns)
    for i = 1:nc,j = 1:ns
        previousGamma = i == 1 ? surface_circulation[j] : wake[i-1,j].gamma
        currentGamma = wake[i,j].gamma
        if symmetric && on_symmetry_plane(top_left(wake[nc,ns]),top_right(wake[nc,ns]))
            gammaH[i,j] = 0.0
        else
            gammaH[i,j] = currentGamma - previousGamma
        end
    end

    # vertical bound vortex circulation strength
    gammaV = Matrix{TF}(undef,nc,ns + 1)
    for i = 1:nc
        currentGamma = wake[i,1].gamma
        if symmetric && on_symmetry_plane(bottom_left(wake[i,1]),top_left(wake[i,1]))
            gammaV[i,1] = 0.0
        else
            gammaV[i,1] = -currentGamma
        end
        for j = 2:ns
            previousGamma = currentGamma
            currentGamma = wake[i,j].gamma
            if symmetric && on_symmetry_plane(bottom_left(wake[i,j]),top_left(wake[i,j]))
                gammaV[i,j] = 0.0
            else
                gammaV[i,j] = previousGamma - currentGamma
            end
        end
        previousGamma = currentGamma
        if symmetric && on_symmetry_plane(bottom_right(wake[i,end]),top_right(wake[i,end]))
            gammaV[i,end] = 0.0
        else
            gammaV[i,end] = previousGamma
        end
    end

    # horizontal bound vortices
    vtk_grid(vtmfile,points,linesH) do vtkfile

        # add metadata
        for (key,value) in pairs(metadata)
            vtkfile[string(key)] = value
        end

        # circulation strength
        vtkfile["circulation"] = reshape(gammaH,:)
    end

    # vertical bound vortices
    vtk_grid(vtmfile,points,linesV) do vtkfile

        # add metadata
        for (key,value) in pairs(metadata)
            vtkfile[string(key)] = value
        end

        # circulation strength
        vtkfile["circulation"] = reshape(gammaV,:)
    end

    if trailing_vortices

        # trailing vortices points as a grid
        xyzT = Array{TF}(undef,3,2,ns + 1)
        for j = 1:ns
            xyzT[:,1,j] = bottom_left(wake[end,j])
            xyzT[:,2,j] = xyzT[:,1,j] + wake_length * xhat
        end
        xyzT[:,1,end] = bottom_right(wake[end,end])
        xyzT[:,2,end] = xyzT[:,1,end] + wake_length * xhat

        # generate points
        pointsT = reshape(xyzT,3,:)

        li = LinearIndices((2,ns + 1))

        # horizontal bound vortices
        linesH = [MeshCell(PolyData.Lines(),[li[1,j],li[1,j+1]]) for j = 1:ns]

        # vertical bound vortices
        linesV = [MeshCell(PolyData.Lines(),[li[1,j],li[2,j]]) for j = 1:ns+1]

        # combine
        linesT = vcat(linesH,linesV)

        # horizontal bound vortex circulation strength
        gammaH = zeros(ns)

        # vertical bound vortex circulation strength
        gammaV = Vector{TF}(undef,ns + 1)
        currentGamma = wake[end,1].gamma
        if symmetric && on_symmetry_plane(bottom_left(wake[end,1]))
            gammaV[1] = 0.0
        else
            gammaV[1] = -currentGamma
        end
        for j = 2:ns
            previousGamma = currentGamma
            currentGamma = wake[end,j].gamma
            if symmetric && on_symmetry_plane(bottom_left(wake[end,j]))
                gammaV[j] = 0.0
            else
                gammaV[j] = previousGamma - currentGamma
            end
        end
        previousGamma = currentGamma
        if symmetric && on_symmetry_plane(bottom_right(wake[end,end]))
            gammaV[j] = 0.0
        else
            gammaV[end] = previousGamma
        end

        # combine
        gammaT = vcat(gammaH,gammaV)

    else

        # generate points
        pointsT = Matrix{TF}(undef,3,ns + 1)
        for j = 1:ns
            pointsT[:,j] = bottom_left(wake[end,j])
        end
        pointsT[:,end] = bottom_right(wake[end,end])

        # horizontal bound vortices
        linesT = [MeshCell(PolyData.Lines(),j:j+1) for j = 1:ns]

        # horizontal bound vortex circulation strength
        gammaT = Vector{TF}(undef,ns)
        for j = 1:ns
            if symmetric && on_symmetry_plane(bottom_left(wake[end,j]),bottom_right(wake[end,j]))
                gammaT[j] = 0.0
            else
                gammaT[j] = -wake[end,j].gamma
            end
        end

    end

    # trailing edge and trailing vortices
    vtk_grid(vtmfile,pointsT,linesT) do vtkfile

        # add metadata
        for (key,value) in pairs(metadata)
            vtkfile[string(key)] = value
        end

        # circulation strength
        vtkfile["circulation"] = reshape(gammaT,:)
    end

    return nothing
end

