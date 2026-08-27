# Legacy near-field VortexLattice force decomposition retained for reference while the
# stability-derivative port is completed. Production near-field loads are
# implemented in `nearfield.jl` and dispatched through
# `near_field_forces!`.

"""
    legacy_near_field_forces!(properties, surfaces, wakes, reference, freestream, Γ;
        dΓdt, additional_velocity, Vh, Vv, symmetric, nwake, surface_id,
        wake_finite_core, wake_shedding_locations, trailing_vortices, xhat)

Calculate local panel forces in the body frame.
"""
function legacy_near_field_forces!(props, surfaces, wakes, ref, fs, Γ;
    dΓdt, additional_velocity, Vh, Vv, symmetric, nwake, surface_id,
    wake_finite_core, wake_shedding_locations, trailing_vortices, xhat,
    interaction_id = surface_id,
    interaction::Bool = true)

    nsurf = length(surfaces)
    TF = eltype(Γ)
    
    # Initialize vectors to hold forces for each surface
    chord_seg_forces = Vector{Matrix{SVector{3, TF}}}()
    span_seg_forces = Vector{Matrix{SVector{3, TF}}}()
    unsteady_forces = Vector{Matrix{SVector{3, TF}}}()

    # loop through receiving surfaces
    iΓ = 0 # index for accessing Γ
    for isurf = 1:nsurf

        receiving = surfaces[isurf]
        nr = length(receiving)
        nr1, nr2 = size(receiving)
        cr = CartesianIndices(receiving)

        # Initialize force matrices for current surface
        current_chord = fill(zero(SVector{3, TF}), nr1, nr2 + 1)
        current_span = fill(zero(SVector{3, TF}), nr1 + 1, nr2)
        current_unsteady = fill(zero(SVector{3, TF}), nr1, nr2)

        # loop through receiving panels
        for i = 1:length(receiving)

            # get panel cartesian index
            I = cr[i]

            # --- Calculate forces on the panel bound vortex --- #

            # bound vortex location
            rc = top_center(receiving[I])

            # freestream velocity
            Vi = freestream_velocity(fs)

            # rotational velocity
            Vi += rotational_velocity(rc, fs, ref)

            # additional velocity field
            if !isnothing(additional_velocity)
                Vi += additional_velocity(rc)
            end

            # velocity due to surface motion
            if !isnothing(Vh)
                Vi += Vh[isurf][i]
            end
            V_streamwise = deepcopy(Vi)

            # induced velocity from surfaces and wakes
            jΓ = 0 # index for accessing Γ
            for jsurf = 1:nsurf

                # number of panels on sending surface
                sending = surfaces[jsurf]
                Ns = length(sending)
                same_interaction_group = interaction || (interaction_id[isurf] == interaction_id[jsurf])

                # see if wake panels are being used
                wake_panels = nwake[jsurf] > 0

                # check if we need to shift shedding locations
                if isnothing(wake_shedding_locations)
                    shedding_locations = nothing
                else
                    shedding_locations = wake_shedding_locations[jsurf]
                end

                same_interaction_group = interaction || (interaction_id[isurf] == interaction_id[jsurf])
                if !same_interaction_group
                    jΓ += Ns
                    continue
                end

                # extract circulation values corresonding to the sending surface
                vΓ = view(Γ, jΓ+1:jΓ+Ns)

                if same_interaction_group
                    # induced velocity from this surface
                    if isurf == jsurf
                    # induced velocity on self
                    Vi += induced_velocity(I, surfaces[jsurf], vΓ;
                        finite_core = surface_id[isurf] != surface_id[jsurf],
                        wake_shedding_locations = shedding_locations,
                        symmetric = symmetric[jsurf],
                        trailing_vortices = trailing_vortices[jsurf] && !wake_panels,
                        xhat = xhat)

                    # streamwise velocity
                    V_streamwise += induced_velocity(I, surfaces[jsurf], vΓ;
                        finite_core = surface_id[isurf] != surface_id[jsurf],
                        wake_shedding_locations = shedding_locations,
                        symmetric = symmetric[jsurf],
                        trailing_vortices = trailing_vortices[jsurf] && !wake_panels,
                        xhat = xhat, skip_leading_edge = true, skip_inside_edges = true, skip_trailing_edge = true)
                else
                    # induced velocity on another surface
                    Vi += induced_velocity(rc, surfaces[jsurf], vΓ;
                        finite_core = surface_id[isurf] != surface_id[jsurf],
                        wake_shedding_locations = shedding_locations,
                        symmetric = symmetric[jsurf],
                        trailing_vortices = trailing_vortices[jsurf] && !wake_panels,
                        xhat = xhat)

                    V_streamwise += induced_velocity(rc, surfaces[jsurf], vΓ;
                        finite_core = surface_id[isurf] != surface_id[jsurf],
                        wake_shedding_locations = shedding_locations,
                        symmetric = symmetric[jsurf],
                        trailing_vortices = trailing_vortices[jsurf] && !wake_panels,
                        xhat = xhat, skip_leading_edge = true, skip_inside_edges = true, skip_trailing_edge = true)
                    end
                end

                # induced velocity from corresponding wake
                if same_interaction_group && wake_panels
                    Vi += induced_velocity(rc, wakes[jsurf];
                        finite_core = wake_finite_core[jsurf] || (surface_id[isurf] != surface_id[jsurf]),
                        symmetric = symmetric[jsurf],
                        nc = nwake[jsurf],
                        trailing_vortices = trailing_vortices[jsurf],
                        xhat = xhat)
#
                    V_streamwise += induced_velocity(rc, wakes[jsurf];
                        finite_core = wake_finite_core[jsurf] || (surface_id[isurf] != surface_id[jsurf]),
                        symmetric = symmetric[jsurf],
                        nc = nwake[jsurf],
                        trailing_vortices = trailing_vortices[jsurf],
                        xhat = xhat, skip_leading_edge = true, skip_inside_edges = true, skip_trailing_edge = true)
                end

                    jΓ += Ns
            end

            # steady part of Kutta-Joukowski theorem
            Γi = I[1] == 1 ? Γ[iΓ+i] : Γ[iΓ+i] - Γ[iΓ+i-1] # net circulation
            Δs = top_vector(receiving[I]) # bound vortex vector
            #Δs = top_center(receiving[I]) - controlpoint(receiving[I]) # bound vortex vector
            tmp = cross(Vi, Δs)
            Fbi = ref.rho*Γi*tmp

            if !isnothing(dΓdt)
                # unsteady part of Kutta-Joukowski theorem

                #TODO: decide whether to divide by perpindicular velocity like
                # Drela does in ASWING?

                dΓdti = I[1] == 1 ? dΓdt[iΓ+i] : (dΓdt[iΓ+i] + dΓdt[iΓ+i-1])/2
                #dΓdti = dΓdt[iΓ+i]
                c = receiving[I].chord
                Fbi += ref.rho*dΓdti*c*tmp
                #ncp = normal(receiving[I])
                #ΔS = (c/nc)*(ref.b/ns)
                #Fbi += ref.rho*dΓdti*(ΔS)*ncp

            end

            # --- Calculate forces on the left bound vortex --- #

            # bound vortex location
            rc = left_center(receiving[I])

            # freestream velocity
            Veff = freestream_velocity(fs)

            # rotational velocity
            Veff += rotational_velocity(rc, fs, ref)

            # additional velocity field
            if !isnothing(additional_velocity)
                Veff += additional_velocity(rc)
            end

            # velocity due to surface motion
            if !isnothing(Vv)
                Veff += Vv[isurf][I[1], I[2]]
            end

            # NOTE: We don't include induced velocity in the effective velocity
            # for the vertical segments because its influence is likely negligible
            # once we take the cross product with the bound vortex vector. This
            # is also assumed in AVL. This could change in the future.

            # steady part of Kutta-Joukowski theorem
            #Γli = I[2] == 1 ? Γ[iΓ+i] : Γ[iΓ+i] - Γ[iΓ+i-1]
            Γli = Γ[iΓ+i]
            #Γi = I[1] == 1 ? Γ[iΓ+i] : Γ[iΓ+i] - Γ[iΓ+i-1] # net circulation
            Δs = left_vector(receiving[I])
            #Δs = controlpoint(receiving[I]) - left_center(receiving[I])
            Fbli = ref.rho*Γli*cross(Veff, Δs)

            # --- Calculate forces on the right bound vortex --- #

            rc = right_center(receiving[I])

            # freestream velocity
            Veff = freestream_velocity(fs)

            # rotational velocity
            Veff += rotational_velocity(rc, fs, ref)

            # additional velocity field
            if !isnothing(additional_velocity)
                Veff += additional_velocity(rc)
            end

            # velocity due to surface motion
            if !isnothing(Vv)
                Veff += Vv[isurf][I[1], I[2]+1]
            end

            # NOTE: We don't include induced velocity in the effective velocity
            # for the vertical segments because its influence is likely negligible
            # once we take the cross product with the bound vortex vector. This
            # is also assumed in AVL. This could change in the future.

            # steady part of Kutta-Joukowski theorem
            Γri = Γ[iΓ+i]
            #Γri = I[2] == 20 ? Γ[iΓ+i] : Γ[iΓ+i] - Γ[iΓ+i-1]
            #Γi = I[1] == 1 ? Γ[iΓ+i] : Γ[iΓ+i] - Γ[iΓ+i-1] # net circulation
            Δs = right_vector(receiving[I])
            #Δs = right_center(receiving[I]) - controlpoint(receiving[I])
            Fbri = ref.rho*Γri*cross(Veff, Δs)

            # store panel circulation, velocity, and forces
            q = 1/2*ref.rho*ref.V^2

            #props[isurf][i] = PanelProperties(Γ[iΓ+i]/ref.V, Vi/ref.V,
            #    Fbi, Fbli, Fbri, V_streamwise)
            
            props[isurf][i] = PanelProperties(Γ[iΓ+i]/ref.V, Vi/ref.V,
                Fbi/(q*ref.S), Fbli/(q*ref.S), Fbri/(q*ref.S), V_streamwise)


            # FORCE CALCULATIONS FOLLOWING IMPERIAL COLLEGE C++ UVLM CODE
            # spanwise segments
            rc = top_center(receiving[I])

            # freestream velocity
            Vi = freestream_velocity(fs)

            # rotational velocity
            Vi += rotational_velocity(rc, fs, ref)

            # additional velocity field
            if !isnothing(additional_velocity)
                Vi += additional_velocity(rc)
            end

            # velocity due to surface motion
            if !isnothing(Vh)
                Vi += Vh[isurf][i]
            end
            V_streamwise = deepcopy(Vi)

            # induced velocity from surfaces and wakes
            jΓ = 0 # index for accessing Γ
            for jsurf = 1:nsurf
#
                # number of panels on sending surface
                sending = surfaces[jsurf]
                Ns = length(sending)
#
                # see if wake panels are being used
                wake_panels = nwake[jsurf] > 0
#
                # check if we need to shift shedding locations
                if isnothing(wake_shedding_locations)
                #    shedding_locations = nothing
                else
                    shedding_locations = wake_shedding_locations[jsurf]
                end
#
                # extract circulation values corresonding to the sending surface
                vΓ = view(Γ, jΓ+1:jΓ+Ns)
#
                same_interaction_group = interaction || (interaction_id[isurf] == interaction_id[jsurf])
                if !same_interaction_group
                    jΓ += Ns
                    continue
                end

                # induced velocity from this surface
                if isurf == jsurf
                    # induced velocity on self
                    Vi += induced_velocity(I, surfaces[jsurf], vΓ;
                        finite_core = surface_id[isurf] != surface_id[jsurf],
                        wake_shedding_locations = shedding_locations,
                        symmetric = symmetric[jsurf],
                        trailing_vortices = trailing_vortices[jsurf] && !wake_panels,
                        xhat = xhat)
#
                    # streamwise velocity
                    V_streamwise += induced_velocity(I, surfaces[jsurf], vΓ;
                        finite_core = surface_id[isurf] != surface_id[jsurf],
                        wake_shedding_locations = shedding_locations,
                        symmetric = symmetric[jsurf],
                        trailing_vortices = trailing_vortices[jsurf] && !wake_panels,
                        xhat = xhat, skip_leading_edge = true, skip_inside_edges = true, skip_trailing_edge = true)
                else
                    # induced velocity on another surface
                    Vi += induced_velocity(rc, surfaces[jsurf], vΓ;
                        finite_core = surface_id[isurf] != surface_id[jsurf],
                        wake_shedding_locations = shedding_locations,
                        symmetric = symmetric[jsurf],
                        trailing_vortices = trailing_vortices[jsurf] && !wake_panels,
                        xhat = xhat)
#
                    V_streamwise += induced_velocity(rc, surfaces[jsurf], vΓ;
                        finite_core = surface_id[isurf] != surface_id[jsurf],
                        wake_shedding_locations = shedding_locations,
                        symmetric = symmetric[jsurf],
                        trailing_vortices = trailing_vortices[jsurf] && !wake_panels,
                        xhat = xhat, skip_leading_edge = true, skip_inside_edges = true, skip_trailing_edge = true)
                end
#
                # induced velocity from corresponding wake
                if same_interaction_group && wake_panels
                    Vi += induced_velocity(rc, wakes[jsurf];
                        finite_core = wake_finite_core[jsurf] || (surface_id[isurf] != surface_id[jsurf]),
                        symmetric = symmetric[jsurf],
                        nc = nwake[jsurf],
                        trailing_vortices = trailing_vortices[jsurf],
                        xhat = xhat)
#
                    V_streamwise += induced_velocity(rc, wakes[jsurf];
                        finite_core = wake_finite_core[jsurf] || (surface_id[isurf] != surface_id[jsurf]),
                        symmetric = symmetric[jsurf],
                        nc = nwake[jsurf],
                        trailing_vortices = trailing_vortices[jsurf],
                        xhat = xhat, skip_leading_edge = true, skip_inside_edges = true, skip_trailing_edge = true)
                end
#
                    jΓ += Ns
            end

            dl = top_vector(receiving[I])
            
            if I[1] == 1
                delta_gamma = Γ[iΓ + i];#-Γ[iΓ + i];
            #} else if (i_M == M){
            #    // Might be needed if TE forces are computed
            #    delta_gamma = lifting_surfaces.gamma[i_surf](i_M-1, i_N);
            else 
                delta_gamma = -Γ[iΓ + i - 1] + Γ[iΓ + i];#Γ[iΓ + i - 1] - Γ[iΓ + i];
            end

            f = ref.rho*delta_gamma*cross(Vi,dl)

            current_span[I[1], I[2]] = f
                
            # trailing edge forces considered zero here
            if I[1] == nr1
                rc = bottom_center(receiving[I])

                # freestream velocity
                Vi = freestream_velocity(fs)
    
                # rotational velocity
                Vi += rotational_velocity(rc, fs, ref)
    
                # additional velocity field
                if !isnothing(additional_velocity)
                    Vi += additional_velocity(rc)
                end
    
                # velocity due to surface motion
                if !isnothing(Vh)
                    Vi += Vh[isurf][i]
                end
                V_streamwise = deepcopy(Vi)
    
                # induced velocity from surfaces and wakes
                jΓ = 0 # index for accessing Γ
                for jsurf = 1:nsurf
    #
                    # number of panels on sending surface
                    sending = surfaces[jsurf]
                    Ns = length(sending)
    #
                    # see if wake panels are being used
                    wake_panels = nwake[jsurf] > 0
    #
                    # check if we need to shift shedding locations
                    if isnothing(wake_shedding_locations)
                        shedding_locations = nothing
                    else
                        shedding_locations = wake_shedding_locations[jsurf]
                    end
    #
                    # extract circulation values corresonding to the sending surface
                    vΓ = view(Γ, jΓ+1:jΓ+Ns)
    #
                    same_interaction_group = interaction || (interaction_id[isurf] == interaction_id[jsurf])
                    if !same_interaction_group
                        jΓ += Ns
                        continue
                    end

                    # induced velocity from this surface
                    if isurf == jsurf
                        # induced velocity on self
                        Vi += induced_velocity(I, surfaces[jsurf], vΓ;
                            finite_core = surface_id[isurf] != surface_id[jsurf],
                            wake_shedding_locations = shedding_locations,
                            symmetric = symmetric[jsurf],
                            trailing_vortices = trailing_vortices[jsurf] && !wake_panels,
                            xhat = xhat)
    #
                        # streamwise velocity
                        V_streamwise += induced_velocity(I, surfaces[jsurf], vΓ;
                            finite_core = surface_id[isurf] != surface_id[jsurf],
                            wake_shedding_locations = shedding_locations,
                            symmetric = symmetric[jsurf],
                            trailing_vortices = trailing_vortices[jsurf] && !wake_panels,
                            xhat = xhat, skip_leading_edge = true, skip_inside_edges = true, skip_trailing_edge = true)
                    else
                        # induced velocity on another surface
                        Vi += induced_velocity(rc, surfaces[jsurf], vΓ;
                            finite_core = surface_id[isurf] != surface_id[jsurf],
                            wake_shedding_locations = shedding_locations,
                            symmetric = symmetric[jsurf],
                            trailing_vortices = trailing_vortices[jsurf] && !wake_panels,
                            xhat = xhat)
    #
                        V_streamwise += induced_velocity(rc, surfaces[jsurf], vΓ;
                            finite_core = surface_id[isurf] != surface_id[jsurf],
                            wake_shedding_locations = shedding_locations,
                            symmetric = symmetric[jsurf],
                            trailing_vortices = trailing_vortices[jsurf] && !wake_panels,
                            xhat = xhat, skip_leading_edge = true, skip_inside_edges = true, skip_trailing_edge = true)
                    end
    #
                    # induced velocity from corresponding wake
                    if same_interaction_group && wake_panels
                        Vi += induced_velocity(rc, wakes[jsurf];
                            finite_core = wake_finite_core[jsurf] || (surface_id[isurf] != surface_id[jsurf]),
                            symmetric = symmetric[jsurf],
                            nc = nwake[jsurf],
                            trailing_vortices = trailing_vortices[jsurf],
                            xhat = xhat)
#
                        V_streamwise += induced_velocity(rc, wakes[jsurf];
                            finite_core = wake_finite_core[jsurf] || (surface_id[isurf] != surface_id[jsurf]),
                            symmetric = symmetric[jsurf],
                            nc = nwake[jsurf],
                            trailing_vortices = trailing_vortices[jsurf],
                            xhat = xhat, skip_leading_edge = true, skip_inside_edges = true, skip_trailing_edge = true)
                    end
    #
                    jΓ += Ns
                end


                dl = bottom_vector(receiving[I])
                delta_gamma = -Γ[iΓ + i]#-Γ[iΓ + i]
                f = ref.rho*delta_gamma*cross(Vi,dl)
                current_span[I[1] + 1, I[2]] = [0.0,0.0,0.0]
            end

            # chord segement forces
            rc = left_center(receiving[I])

            # freestream velocity
            Vi = freestream_velocity(fs)

            # rotational velocity
            Vi += rotational_velocity(rc, fs, ref)

            # additional velocity field
            if !isnothing(additional_velocity)
                Vi += additional_velocity(rc)
            end

            # velocity due to surface motion
            if !isnothing(Vv)
                Vi += Vv[isurf][I[1], I[2]]
            end
            V_streamwise = deepcopy(Vi)

            # induced velocity from surfaces and wakes
            jΓ = 0 # index for accessing Γ
            for jsurf = 1:nsurf
###
                # number of panels on sending surface
                sending = surfaces[jsurf]
                Ns = length(sending)
###
                # see if wake panels are being used
                wake_panels = nwake[jsurf] > 0
###
                # check if we need to shift shedding locations
                if isnothing(wake_shedding_locations)
                    shedding_locations = nothing
                else
                    shedding_locations = wake_shedding_locations[jsurf]
                end
###
                # extract circulation values corresonding to the sending surface
                vΓ = view(Γ, jΓ+1:jΓ+Ns)
###
                same_interaction_group = interaction || (interaction_id[isurf] == interaction_id[jsurf])
                if !same_interaction_group
                    jΓ += Ns
                    continue
                end
#
                # induced velocity from this surface
                if isurf == jsurf
                    if wake_panels
                        Vi += induced_velocity(rc, wakes[jsurf];
                            finite_core = wake_finite_core[jsurf],
                            symmetric = symmetric[jsurf],
                            nc = nwake[jsurf],
                            trailing_vortices = trailing_vortices[jsurf],
                            xhat = xhat)
#
                        V_streamwise += induced_velocity(rc, wakes[jsurf];
                            finite_core = wake_finite_core[jsurf],
                            symmetric = symmetric[jsurf],
                            nc = nwake[jsurf],
                            trailing_vortices = trailing_vortices[jsurf],
                            xhat = xhat, skip_leading_edge = true, skip_inside_edges = true, skip_trailing_edge = true)
                    end
                    jΓ += Ns
                    continue
                else
                    # induced velocity on another surface
                    Vi += induced_velocity(rc, surfaces[jsurf], vΓ;
                        finite_core = surface_id[isurf] != surface_id[jsurf],
                        wake_shedding_locations = shedding_locations,
                        symmetric = symmetric[jsurf],
                        trailing_vortices = trailing_vortices[jsurf] && !wake_panels,
                        xhat = xhat)
###
                    V_streamwise += induced_velocity(rc, surfaces[jsurf], vΓ;
                        finite_core = surface_id[isurf] != surface_id[jsurf],
                        wake_shedding_locations = shedding_locations,
                        symmetric = symmetric[jsurf],
                        trailing_vortices = trailing_vortices[jsurf] && !wake_panels,
                        xhat = xhat, skip_leading_edge = true, skip_inside_edges = true, skip_trailing_edge = true)
                end
###
                # induced velocity from corresponding wake
                if same_interaction_group && wake_panels
                    Vi += induced_velocity(rc, wakes[jsurf];
                        finite_core = wake_finite_core[jsurf] || (surface_id[isurf] != surface_id[jsurf]),
                        symmetric = symmetric[jsurf],
                        nc = nwake[jsurf],
                        trailing_vortices = trailing_vortices[jsurf],
                        xhat = xhat)
##
                    V_streamwise += induced_velocity(rc, wakes[jsurf];
                        finite_core = wake_finite_core[jsurf] || (surface_id[isurf] != surface_id[jsurf]),
                        symmetric = symmetric[jsurf],
                        nc = nwake[jsurf],
                        trailing_vortices = trailing_vortices[jsurf],
                        xhat = xhat, skip_leading_edge = true, skip_inside_edges = true, skip_trailing_edge = true)
                end
###
                    jΓ += Ns
            end

            dl = left_vector(receiving[I]);

            if I[2] == 1
                delta_gamma = Γ[iΓ + i];
            else 
                delta_gamma = Γ[iΓ + i] - Γ[iΓ + i - nr1];
            end
            

            f = ref.rho*delta_gamma*cross(Vi, dl);
            current_chord[I[1], I[2]] = f
            
            # Consider last collum of chordwise panels
            # Final chordwise column
        if I[2] == nr2
            rc = right_center(receiving[I])

                        # freestream velocity
            Vi = freestream_velocity(fs)

            # rotational velocity
            Vi += rotational_velocity(rc, fs, ref)

            # additional velocity field
            if !isnothing(additional_velocity)
                Vi += additional_velocity(rc)
            end

            # velocity due to surface motion
            if !isnothing(Vv)
                Vi += Vv[isurf][I[1], I[2]+1]
            end
            V_streamwise = deepcopy(Vi)

            # induced velocity from surfaces and wakes
            jΓ = 0 # index for accessing Γ
            for jsurf = 1:nsurf
##
                # number of panels on sending surface
                sending = surfaces[jsurf]
                Ns = length(sending)
##
                # see if wake panels are being used
                wake_panels = nwake[jsurf] > 0
##
                # check if we need to shift shedding locations
                if isnothing(wake_shedding_locations)
                    shedding_locations = nothing
                else
                    shedding_locations = wake_shedding_locations[jsurf]
                end
##
                # extract circulation values corresonding to the sending surface
                vΓ = view(Γ, jΓ+1:jΓ+Ns)
##
                same_interaction_group = interaction || (interaction_id[isurf] == interaction_id[jsurf])
                if !same_interaction_group
                    jΓ += Ns
                    continue
                end
#
                # induced velocity from this surface
                if isurf == jsurf
                    if wake_panels
                        Vi += induced_velocity(rc, wakes[jsurf];
                            finite_core = wake_finite_core[jsurf],
                            symmetric = symmetric[jsurf],
                            nc = nwake[jsurf],
                            trailing_vortices = trailing_vortices[jsurf],
                            xhat = xhat)
#
                        V_streamwise += induced_velocity(rc, wakes[jsurf];
                            finite_core = wake_finite_core[jsurf],
                            symmetric = symmetric[jsurf],
                            nc = nwake[jsurf],
                            trailing_vortices = trailing_vortices[jsurf],
                            xhat = xhat, skip_leading_edge = true, skip_inside_edges = true, skip_trailing_edge = true)
                    end
                    jΓ += Ns
                    continue
                else
                    # induced velocity on another surface
                    Vi += induced_velocity(rc, surfaces[jsurf], vΓ;
                        finite_core = surface_id[isurf] != surface_id[jsurf],
                        wake_shedding_locations = shedding_locations,
                        symmetric = symmetric[jsurf],
                        trailing_vortices = trailing_vortices[jsurf] && !wake_panels,
                        xhat = xhat)
##
                    V_streamwise += induced_velocity(rc, surfaces[jsurf], vΓ;
                        finite_core = surface_id[isurf] != surface_id[jsurf],
                        wake_shedding_locations = shedding_locations,
                        symmetric = symmetric[jsurf],
                        trailing_vortices = trailing_vortices[jsurf] && !wake_panels,
                        xhat = xhat, skip_leading_edge = true, skip_inside_edges = true, skip_trailing_edge = true)
                end
##
                # induced velocity from corresponding wake
                if same_interaction_group && wake_panels
                    Vi += induced_velocity(rc, wakes[jsurf];
                        finite_core = wake_finite_core[jsurf] || (surface_id[isurf] != surface_id[jsurf]),
                        symmetric = symmetric[jsurf],
                        nc = nwake[jsurf],
                        trailing_vortices = trailing_vortices[jsurf],
                        xhat = xhat)
##
                    V_streamwise += induced_velocity(rc, wakes[jsurf];
                        finite_core = wake_finite_core[jsurf] || (surface_id[isurf] != surface_id[jsurf]),
                        symmetric = symmetric[jsurf],
                        nc = nwake[jsurf],
                        trailing_vortices = trailing_vortices[jsurf],
                        xhat = xhat, skip_leading_edge = true, skip_inside_edges = true, skip_trailing_edge = true)
                end
##
                    jΓ += Ns
            end
            dl = right_vector(receiving[I])

            delta_gamma = Γ[iΓ + i]#-Γ[iΓ + i]
            f = ref.rho * delta_gamma * cross(Vi, dl)
            #for comp in 1:3
                current_chord[I[1], I[2] + 1] = f
            #end


        end
        if !isnothing(dΓdt)
            # unsteady loads
            normi = normal(receiving[I])
            Δs_span = norm(top_vector(receiving[I]))
            Δs_chord = norm(left_vector(receiving[I]))
            area = Δs_span * Δs_chord
            current_unsteady[I[1],I[2]] = ref.rho*area*normi*dΓdt[iΓ+i]
        end
    end


        # Store forces for current surface
        push!(chord_seg_forces, current_chord)
        push!(span_seg_forces, current_span)
        push!(unsteady_forces, current_unsteady)

        # increment Γ index for receiving panels
        iΓ += nr
    end

    return props, chord_seg_forces, span_seg_forces, unsteady_forces
end


## ==============================================================================
## REF. NEARFIELD.JL - CORRECTED (NO SINGULARITY)
## ==============================================================================
#
#function near_field_forces!(props, surfaces, wakes, ref, fs, Γ;
#    dΓdt, additional_velocity, Vh, Vv, symmetric, nwake, surface_id,
#    wake_finite_core, wake_shedding_locations, trailing_vortices, xhat)
#
#    nsurf = length(surfaces)
#    TF = eltype(Γ)
#    
#    # Initialize output containers
#    chord_seg_forces = Vector{Matrix{SVector{3, TF}}}()
#    span_seg_forces = Vector{Matrix{SVector{3, TF}}}()
#    unsteady_forces = Vector{Matrix{SVector{3, TF}}}()
#
#    iΓ = 0 # index for accessing Γ global vector
#
#    for isurf = 1:nsurf
#        receiving = surfaces[isurf]
#        nr1, nr2 = size(receiving)
#        
#        # Force matrices for current surface
#        current_chord = fill(zero(SVector{3, TF}), nr1, nr2 + 1)
#        current_span = fill(zero(SVector{3, TF}), nr1 + 1, nr2)
#        current_unsteady = fill(zero(SVector{3, TF}), nr1, nr2)
#
#        for i = 1:length(receiving)
#            I = CartesianIndices(receiving)[i] # [chord_idx, span_idx]
#
#            # -------------------------------------------------------
#            # 1. COMPUTE FLUID VELOCITY (V_fluid)
#            # -------------------------------------------------------
#            rc = top_center(receiving[I])
#            
#            # A. Freestream + Rotation + Gust
#            V_fluid = freestream_velocity(fs) + rotational_velocity(rc, fs, ref)
#            if !isnothing(additional_velocity)
#                V_fluid += additional_velocity(rc)
#            end
#            
#            # B. Induced Velocity
#            jΓ = 0
#            for jsurf = 1:nsurf
#                sending = surfaces[jsurf]
#                Ns = length(sending)
#                vΓ = view(Γ, jΓ+1:jΓ+Ns)
#                
#                use_finite_core = surface_id[isurf] != surface_id[jsurf]
#                
#                # Induced by Surface
#                # CRITICAL FIX: Use 'I' when isurf==jsurf to skip self-induction
#                if isurf == jsurf
#                    V_fluid += induced_velocity(I, surfaces[jsurf], vΓ;
#                        finite_core = use_finite_core,
#                        wake_shedding_locations = isnothing(wake_shedding_locations) ? nothing : wake_shedding_locations[jsurf],
#                        symmetric = symmetric[jsurf],
#                        trailing_vortices = trailing_vortices[jsurf] && (nwake[jsurf] == 0),
#                        xhat = xhat)
#                else
#                    V_fluid += induced_velocity(rc, surfaces[jsurf], vΓ;
#                        finite_core = use_finite_core,
#                        wake_shedding_locations = isnothing(wake_shedding_locations) ? nothing : wake_shedding_locations[jsurf],
#                        symmetric = symmetric[jsurf],
#                        trailing_vortices = trailing_vortices[jsurf] && (nwake[jsurf] == 0),
#                        xhat = xhat)
#                end
#                
#                # Induced by Wake
#                if nwake[jsurf] > 0
#                    V_fluid += induced_velocity(rc, wakes[jsurf];
#                        finite_core = wake_finite_core[jsurf] || use_finite_core,
#                        symmetric = symmetric[jsurf],
#                        nc = nwake[jsurf],
#                        trailing_vortices = trailing_vortices[jsurf],
#                        xhat = xhat)
#                end
#                jΓ += Ns
#            end
#            
#            # -------------------------------------------------------
#            # 2. RELATIVE VELOCITY (V_rel = V_fluid - V_grid)
#            # -------------------------------------------------------
#            # Vh contains (r_old - r_new)/dt, which is -V_grid. So we ADD it.
#            V_grid_neg = isnothing(Vh) ? zero(SVector{3, TF}) : Vh[isurf][i]
#            V_rel = V_fluid + V_grid_neg 
#            V_streamwise = deepcopy(V_fluid) 
#
#            # -------------------------------------------------------
#            # 3. KUTTA-JOUKOWSKI FORCES
#            # -------------------------------------------------------
#            
#            # A. Spanwise Force
#            dl_span = top_vector(receiving[I])
#            Γ_local = Γ[iΓ+i]
#            delta_Gamma_span = (I[1] == 1) ? Γ_local : (Γ_local - Γ[iΓ + i - 1])
#            
#            f_span = ref.rho * delta_Gamma_span * cross(V_rel, dl_span)
#            current_span[I[1], I[2]] = f_span
#
#            # B. Chordwise Forces
#            dl_left = left_vector(receiving[I])
#            delta_Gamma_left = (I[2] == 1) ? Γ_local : (Γ_local - Γ[iΓ + i - nr1])
#            f_chord_left = ref.rho * delta_Gamma_left * cross(V_rel, dl_left)
#            current_chord[I[1], I[2]] = f_chord_left
#            
#            if I[2] == nr2
#                dl_right = right_vector(receiving[I])
#                delta_Gamma_right = -Γ_local 
#                f_chord_right = ref.rho * delta_Gamma_right * cross(V_rel, dl_right)
#                current_chord[I[1], I[2] + 1] = f_chord_right
#            end
#            
#            # -------------------------------------------------------
#            # 4. UNSTEADY APPARENT MASS FORCE
#            # -------------------------------------------------------
#            if !isnothing(dΓdt)
#                n_vec = normal(receiving[I])
#                area = receiving[I].chord * norm(dl_span)
#                f_unsteady = ref.rho * area * dΓdt[iΓ+i] * n_vec
#                current_unsteady[I[1], I[2]] = f_unsteady
#            else
#                f_unsteady = zero(SVector{3, TF})
#            end
#
#            # -------------------------------------------------------
#            # 5. STORE PROPERTIES
#            # -------------------------------------------------------
#            q_dyn = 0.5 * ref.rho * ref.V^2
#            
#            props[isurf][i] = PanelProperties(
#                Γ_local / ref.V,
#                V_rel / ref.V,
#                (f_span + f_unsteady) / (q_dyn * ref.S),
#                f_chord_left / (q_dyn * ref.S),
#                (I[2] == nr2 ? f_chord_right : zero(SVector{3, TF})) / (q_dyn * ref.S),
#                V_streamwise
#            )
#        end
#        
#        push!(chord_seg_forces, current_chord)
#        push!(span_seg_forces, current_span)
#        push!(unsteady_forces, current_unsteady)
#        
#        iΓ += length(receiving)
#    end
#
#    return props, chord_seg_forces, span_seg_forces, unsteady_forces
#end

#"""
#    near_field_forces!(...)
#
#Calculates local panel forces in the body frame with corrected relative velocity.
#"""
#function near_field_forces!(props, surfaces, wakes, ref, fs, Γ;
#                            dΓdt, additional_velocity, Vh, Vv, symmetric, nwake, surface_id,
#                            wake_finite_core, wake_shedding_locations, trailing_vortices, xhat)
#
#    nsurf = length(surfaces)
#    TF = eltype(Γ)
#
#    # Initialize vectors to hold forces for each surface
#    chord_seg_forces = Vector{Matrix{SVector{3, TF}}}()
#    span_seg_forces = Vector{Matrix{SVector{3, TF}}}()
#    unsteady_forces = Vector{Matrix{SVector{3, TF}}}()
#
#    iΓ = 0 # index for accessing Γ
#    for isurf = 1:nsurf
#        receiving = surfaces[isurf]
#        nr = length(receiving)
#        nr1, nr2 = size(receiving)
#        cr = CartesianIndices(receiving)
#
#        # Initialize force matrices for the current surface
#        current_chord = fill(zero(SVector{3, TF}), nr1, nr2 + 1)
#        current_span = fill(zero(SVector{3, TF}), nr1 + 1, nr2)
#        current_unsteady = fill(zero(SVector{3, TF}), nr1, nr2)
#
#        for i = 1:length(receiving)
#            I = cr[i] # Cartesian index of the panel
#
#            # --- Calculate Velocity Components ---
#            # This logic is now separated to clearly define V_fluid and V_segment
#
#            # --- 1. Horizontal (Spanwise) Bound Vortex ---
#            rc_h = top_center(receiving[I])
#
#            # Calculate V_fluid at the horizontal vortex segment
#            V_fluid_h = freestream_velocity(fs)
#            V_fluid_h += rotational_velocity(rc_h, fs, ref) # Note: This is zero in your main.jl but kept for generality
#            if !isnothing(additional_velocity)
#                V_fluid_h += additional_velocity(rc_h)
#            end
#
#            # Add induced velocities from all surfaces and wakes
#            jΓ_h = 0
#            for jsurf = 1:nsurf
#                vΓ_h = view(Γ, jΓ_h+1:jΓ_h+length(surfaces[jsurf]))
#                # Induced from bound vortices
#                V_fluid_h += induced_velocity(isurf == jsurf ? I : rc_h, surfaces[jsurf], vΓ_h; symmetric=symmetric[jsurf], xhat=xhat)
#                # Induced from wakes
#                if nwake[jsurf] > 0
#                    V_fluid_h += induced_velocity(rc_h, wakes[jsurf]; nc=nwake[jsurf], symmetric=symmetric[jsurf], xhat=xhat)
#                end
#                jΓ_h += length(surfaces[jsurf])
#            end
#
#            # Get the velocity of the segment itself (from whirl and spin)
#            V_segment_h = !isnothing(Vh) ? Vh[isurf][I] : zero(SVector{3, TF})
#
#            # --- CORRECTED Relative Velocity for Horizontal Segment ---
#            V_relative_h = V_fluid_h - V_segment_h
#
#            # --- 2. Vertical (Chordwise) Bound Vortices (Left and Right) ---
#            # Left segment
#            rc_v_left = left_center(receiving[I])
#            V_segment_v_left = !isnothing(Vv) ? Vv[isurf][I[1], I[2]] : zero(SVector{3, TF})
#            V_fluid_v_left = freestream_velocity(fs) + rotational_velocity(rc_v_left, fs, ref) # Induced velocities are ignored on vertical segments as per original logic
#            if !isnothing(additional_velocity); V_fluid_v_left += additional_velocity(rc_v_left); end
#            V_relative_v_left = V_fluid_v_left - V_segment_v_left
#
#            # Right segment
#            rc_v_right = right_center(receiving[I])
#            V_segment_v_right = !isnothing(Vv) ? Vv[isurf][I[1], I[2]+1] : zero(SVector{3, TF})
#            V_fluid_v_right = freestream_velocity(fs) + rotational_velocity(rc_v_right, fs, ref)
#            if !isnothing(additional_velocity); V_fluid_v_right += additional_velocity(rc_v_right); end
#            V_relative_v_right = V_fluid_v_right - V_segment_v_right
#
#            # --- Calculate Forces using CORRECTED Relative Velocities ---
#            Γi = Γ[iΓ+i]
#            q = 0.5 * ref.rho * ref.V^2
#
#            # Force on horizontal (spanwise) segment
#            delta_gamma_h = (I[1] == 1) ? Γi : (Γi - Γ[iΓ+i-1])
#            dl_h = top_vector(receiving[I])
#            f_span = ref.rho * delta_gamma_h * cross(V_relative_h, dl_h)
#            current_span[I[1], I[2]] = f_span
#
#            # Force on left vertical (chordwise) segment
#            delta_gamma_v_left = (I[2] == 1) ? Γi : (Γi - Γ[iΓ+i-nr1])
#            dl_v_left = left_vector(receiving[I])
#            f_chord_left = ref.rho * delta_gamma_v_left * cross(V_relative_v_left, dl_v_left)
#            current_chord[I[1], I[2]] = f_chord_left
#
#            # Force on right vertical (chordwise) segment (for the last column)
#            if I[2] == nr2
#                delta_gamma_v_right = -Γi # Shed vortex
#                dl_v_right = right_vector(receiving[I])
#                f_chord_right = ref.rho * delta_gamma_v_right * cross(V_relative_v_right, dl_v_right)
#                current_chord[I[1], I[2]+1] = f_chord_right
#            end
#
#            # Trailing edge force (last row)
#            if I[1] == nr1
#                dl_h_te = bottom_vector(receiving[I])
#                delta_gamma_te = -Γi
#                # Using an average velocity for the TE segment
#                V_relative_te = (V_relative_v_left + V_relative_v_right) / 2 # Approximation
#                f_span_te = ref.rho * delta_gamma_te * cross(V_relative_te, dl_h_te)
#                current_span[I[1]+1, I[2]] = f_span_te
#            end
#
#            # Unsteady Force (dGamma/dt term)
#            if !isnothing(dΓdt)
#                area_approx = norm(top_vector(receiving[I])) * norm(left_vector(receiving[I]))
#                normal_vec = normal(receiving[I])
#                current_unsteady[I[1], I[2]] = ref.rho * area_approx * dΓdt[iΓ+i] * normal_vec
#            end
#
#            # Store simplified panel properties for consistency with the rest of the code
#            # Note: cfb, cfl, cfr are now approximations as forces are distributed differently
#            props[isurf][i] = PanelProperties(Γi / ref.V, V_relative_h / ref.V,
#                                              f_span / (q * ref.S), f_chord_left / (q * ref.S),
#                                              (I[2]==nr2 ? current_chord[I[1],I[2]+1] : zero(SVector{3, TF})) / (q * ref.S),
#                                              V_fluid_h)
#        end
#
#        push!(chord_seg_forces, current_chord)
#        push!(span_seg_forces, current_span)
#        push!(unsteady_forces, current_unsteady)
#        iΓ += nr
#    end
#
#    return props, chord_seg_forces, span_seg_forces, unsteady_forces
#end
#
#"""
#    near_field_forces!(...)
#
#Calculates local panel forces based on the Kutta-Joukowski theorem using the
#true relative velocity, in coherence with standard UVLM implementations like DUST.
#"""
#function near_field_forces!(props, surfaces, wakes, ref, fs, Γ;
#                            dΓdt, additional_velocity, Vh, Vv, symmetric, nwake, surface_id,
#                            wake_finite_core, wake_shedding_locations, trailing_vortices, xhat)
#
#    nsurf = length(surfaces)
#    TF = eltype(Γ)
#
#    chord_seg_forces = Vector{Matrix{SVector{3, TF}}}()
#    span_seg_forces = Vector{Matrix{SVector{3, TF}}}()
#    unsteady_forces = Vector{Matrix{SVector{3, TF}}}()
#
#    iΓ = 0
#    for isurf = 1:nsurf
#        receiving = surfaces[isurf]
#        nr1, nr2 = size(receiving)
#        cr = CartesianIndices(receiving)
#
#        current_chord = fill(zero(SVector{3, TF}), nr1, nr2 + 1)
#        current_span = fill(zero(SVector{3, TF}), nr1 + 1, nr2)
#        current_unsteady = fill(zero(SVector{3, TF}), nr1, nr2)
#
#        for i = 1:length(receiving)
#            I = cr[i]
#
#            # --- Calculate Relative Velocity for Each Vortex Segment ---
#
#            # --- 1. Horizontal (Spanwise) Bound Vortex ---
#            rc_h = top_center(receiving[I])
#
#            # 1a. Calculate V_fluid = V_freestream + V_induced
#            V_fluid_h = freestream_velocity(fs)
#            jΓ_h = 0
#            for jsurf = 1:nsurf
#                vΓ_h = view(Γ, jΓ_h+1:jΓ_h+length(surfaces[jsurf]))
#                V_fluid_h += induced_velocity(isurf == jsurf ? I : rc_h, surfaces[jsurf], vΓ_h; symmetric=symmetric[jsurf], xhat=xhat)
#                if nwake[jsurf] > 0
#                    V_fluid_h += induced_velocity(rc_h, wakes[jsurf]; nc=nwake[jsurf], symmetric=symmetric[jsurf], xhat=xhat)
#                end
#                jΓ_h += length(surfaces[jsurf])
#            end
#
#            # 1b. Calculate V_body = V_spin + V_whirl
#            V_spin_h = !isnothing(additional_velocity) ? additional_velocity(rc_h) : zero(SVector{3, TF})
#            V_whirl_h = !isnothing(Vh) ? Vh[isurf][I] : zero(SVector{3, TF})
#            V_body_h = V_spin_h + V_whirl_h
#
#            # 1c. True Relative Velocity
#            V_relative_h = V_fluid_h - V_body_h
#
#            # --- 2. Vertical (Chordwise) Bound Vortices (Approximation) ---
#            # As per original logic, induced velocity is neglected for these smaller segments.
#            rc_v_left = left_center(receiving[I])
#            V_fluid_v_left = freestream_velocity(fs)
#            V_spin_v_left = !isnothing(additional_velocity) ? additional_velocity(rc_v_left) : zero(SVector{3, TF})
#            V_whirl_v_left = !isnothing(Vv) ? Vv[isurf][I[1], I[2]] : zero(SVector{3, TF})
#            V_relative_v_left = V_fluid_v_left - (V_spin_v_left + V_whirl_v_left)
#
#            rc_v_right = right_center(receiving[I])
#            V_fluid_v_right = freestream_velocity(fs)
#            V_spin_v_right = !isnothing(additional_velocity) ? additional_velocity(rc_v_right) : zero(SVector{3, TF})
#            V_whirl_v_right = !isnothing(Vv) ? Vv[isurf][I[1], I[2]+1] : zero(SVector{3, TF})
#            V_relative_v_right = V_fluid_v_right - (V_spin_v_right + V_whirl_v_right)
#
#            # --- Calculate Forces using Correct Relative Velocities ---
#            Γi = Γ[iΓ+i]
#
#            # Force on horizontal (spanwise) segment
#            delta_gamma_h = (I[1] == 1) ? Γi : (Γi - Γ[iΓ+i-1])
#            dl_h = top_vector(receiving[I])
#            f_span = ref.rho * delta_gamma_h * cross(V_relative_h, dl_h)
#            current_span[I[1], I[2]] = f_span
#
#            # Force on left vertical (chordwise) segment
#            delta_gamma_v_left = (I[2] == 1) ? Γi : (Γi - Γ[iΓ+i-nr1])
#            dl_v_left = left_vector(receiving[I])
#            f_chord_left = ref.rho * delta_gamma_v_left * cross(V_relative_v_left, dl_v_left)
#            current_chord[I[1], I[2]] = f_chord_left
#
#            # Force on right vertical (chordwise) segment (last column)
#            if I[2] == nr2
#                delta_gamma_v_right = -Γi
#                dl_v_right = right_vector(receiving[I])
#                f_chord_right = ref.rho * delta_gamma_v_right * cross(V_relative_v_right, dl_v_right)
#                current_chord[I[1], I[2]+1] = f_chord_right
#            end
#            
#            # Force on trailing edge segment (last row)
#            if I[1] == nr1
#                dl_h_te = bottom_vector(receiving[I])
#                delta_gamma_te = -Γi
#                V_relative_te = (V_relative_v_left + V_relative_v_right) / 2 # Approximation
#                f_span_te = ref.rho * delta_gamma_te * cross(V_relative_te, dl_h_te)
#                current_span[I[1]+1, I[2]] = f_span_te
#            end
#
#            # Unsteady Force (dGamma/dt term)
#            if !isnothing(dΓdt)
#                area_approx = norm(top_vector(receiving[I])) * norm(left_vector(receiving[I]))
#                normal_vec = normal(receiving[I])
#                current_unsteady[I[1],I[2]] = ref.rho * area_approx * dΓdt[iΓ+i] * normal_vec
#            end
#
#            # Store panel properties
#            q = 0.5 * ref.rho * ref.V^2
#            props[isurf][i] = PanelProperties(Γi/ref.V, V_relative_h/ref.V, f_span/(q*ref.S), 
#                                              f_chord_left/(q*ref.S), zero(SVector{3, TF}), V_fluid_h)
#
#        end
#
#        push!(chord_seg_forces, current_chord)
#        push!(span_seg_forces, current_span)
#        push!(unsteady_forces, current_unsteady)
#        iΓ += length(receiving)
#    end
#
#    return props, chord_seg_forces, span_seg_forces, unsteady_forces
#end

"""
    legacy_near_field_forces_derivatives!(properties, dproperties, surfaces, reference,
        freestream, Γ, dΓ; dΓdt, additional_velocity, Vh, Vv, symmetric, nwake,
        surface_id, wake_finite_core, wake_shedding_locations, trailing_vortices, xhat)

Version of [`near_field_forces!`](@ref) that also calculates the derivatives of
the local panel forces with respect to the freestream variables.
"""
legacy_near_field_forces_derivatives!

function legacy_near_field_forces_derivatives!(props, dprops, surfaces, wakes,
    ref, fs, Γ, dΓ; dΓdt, additional_velocity, Vh, Vv, symmetric, nwake,
    surface_id, wake_finite_core, wake_shedding_locations, trailing_vortices, xhat,
    interaction_id = surface_id,
    interaction::Bool = true)

    # unpack derivatives
    props_a, props_b, props_p, props_q, props_r = dprops
    Γ_a, Γ_b, Γ_p, Γ_q, Γ_r = dΓ

    # number of surfaces
    nsurf = length(surfaces)

    # loop through receiving surfaces
    iΓ = 0 # index for accessing Γ
    for isurf = 1:nsurf

        receiving = surfaces[isurf]
        nr = length(receiving)
        nr1, nr2 = size(receiving)
        cr = CartesianIndices(receiving)

        # loop through receiving panels
        for i = 1:length(receiving)

            I = cr[i]

            # --- Calculate forces on the panel bound vortex -- #
            rc = top_center(receiving[I])

            # freestream velocity
            Vi, dVi = freestream_velocity_derivatives(fs)
            Vi_a, Vi_b = dVi

            # rotational velocity
            Vrot, dVrot = rotational_velocity_derivatives(rc, fs, ref)
            Vi += Vrot
            Vi_p, Vi_q, Vi_r = dVrot

            # additional velocity field
            if !isnothing(additional_velocity)
                Vi += additional_velocity(rc)
            end

            # velocity due to surface motion
            if !isnothing(Vh)
                Vi += Vh[isurf][i]
            end

            V_streamwise = deepcopy(Vi)

            # induced velocity from surfaces and wakes
            jΓ = 0 # index for accessing Γ
            for jsurf = 1:nsurf

                # number of panels on sending surface
                sending = surfaces[jsurf]
                Ns = length(sending)

                same_interaction_group = interaction || (interaction_id[isurf] == interaction_id[jsurf])
                if !same_interaction_group
                    jΓ += Ns
                    continue
                end

                # see if wake panels are being used
                wake_panels = nwake[jsurf] > 0

                # check if we need to shift shedding locations
                if isnothing(wake_shedding_locations)
                    shedding_locations = nothing
                else
                    shedding_locations = wake_shedding_locations[jsurf]
                end

                # extract circulation values corresonding to the sending surface
                vΓ = view(Γ, jΓ+1:jΓ+Ns)

                vΓ_a = view(Γ_a, jΓ+1:jΓ+Ns)
                vΓ_b = view(Γ_b, jΓ+1:jΓ+Ns)
                vΓ_p = view(Γ_p, jΓ+1:jΓ+Ns)
                vΓ_q = view(Γ_q, jΓ+1:jΓ+Ns)
                vΓ_r = view(Γ_r, jΓ+1:jΓ+Ns)

                vdΓ = (vΓ_a, vΓ_b, vΓ_p, vΓ_q, vΓ_r)

                # induced velocity from this surface
                if isurf == jsurf
                    # induced velocity on self
                    Vind, dVind = induced_velocity_derivatives(I, surfaces[jsurf], vΓ, vdΓ;
                        finite_core = surface_id[isurf] != surface_id[jsurf],
                        wake_shedding_locations = shedding_locations,
                        symmetric = symmetric[jsurf],
                        trailing_vortices = trailing_vortices[jsurf] && !wake_panels,
                        xhat = xhat)

                    Vind_stream, _ = induced_velocity_derivatives(I, surfaces[jsurf], vΓ, vdΓ;
                        finite_core = surface_id[isurf] != surface_id[jsurf],
                        wake_shedding_locations = shedding_locations,
                        symmetric = symmetric[jsurf],
                        trailing_vortices = trailing_vortices[jsurf] && !wake_panels,
                        xhat = xhat, skip_leading_edge = true, skip_inside_edges = true, skip_trailing_edge = true)
                else
                    # induced velocity on another surface
                    Vind, dVind = induced_velocity_derivatives(rc, surfaces[jsurf], vΓ, vdΓ;
                        finite_core = surface_id[isurf] != surface_id[jsurf],
                        wake_shedding_locations = shedding_locations,
                        symmetric = symmetric[jsurf],
                        trailing_vortices = trailing_vortices[jsurf] && !wake_panels,
                        xhat = xhat)

                    Vind_stream, _ = induced_velocity_derivatives(rc, surfaces[jsurf], vΓ, vdΓ;
                        finite_core = surface_id[isurf] != surface_id[jsurf],
                        wake_shedding_locations = shedding_locations,
                        symmetric = symmetric[jsurf],
                        trailing_vortices = trailing_vortices[jsurf] && !wake_panels,
                        xhat = xhat, skip_leading_edge = true, skip_inside_edges = true, skip_trailing_edge = true)
                end

                Vind_a, Vind_b, Vind_p, Vind_q, Vind_r = dVind

                Vi += Vind
                V_streamwise += Vind_stream

                Vi_a += Vind_a
                Vi_b += Vind_b
                Vi_p += Vind_p
                Vi_q += Vind_q
                Vi_r += Vind_r

                # induced velocity from corresponding wake
                if same_interaction_group && wake_panels
                    Vi += induced_velocity(rc, wakes[jsurf];
                        finite_core = wake_finite_core[jsurf] || surface_id[isurf] != surface_id[jsurf],
                        symmetric = symmetric[jsurf],
                        nc = nwake[jsurf],
                        trailing_vortices = trailing_vortices[jsurf],
                        xhat = xhat)

                    V_streamwise += induced_velocity(rc, wakes[jsurf];
                        finite_core = wake_finite_core[jsurf] || surface_id[isurf] != surface_id[jsurf],
                        symmetric = symmetric[jsurf],
                        nc = nwake[jsurf],
                        trailing_vortices = trailing_vortices[jsurf],
                        xhat = xhat, skip_leading_edge = true, skip_inside_edges = true, skip_trailing_edge = true)
                end

                    jΓ += Ns
            end

            # steady part of Kutta-Joukowski theorem
            if I[1] == 1
                Γi = Γ[iΓ+i]

                Γi_a = Γ_a[iΓ+i]
                Γi_b = Γ_b[iΓ+i]
                Γi_p = Γ_p[iΓ+i]
                Γi_q = Γ_q[iΓ+i]
                Γi_r = Γ_r[iΓ+i]
            else
                Γi = Γ[iΓ+i] - Γ[iΓ+i-1]

                Γi_a = Γ_a[iΓ+i] - Γ_a[iΓ+i-1]
                Γi_b = Γ_b[iΓ+i] - Γ_b[iΓ+i-1]
                Γi_p = Γ_p[iΓ+i] - Γ_p[iΓ+i-1]
                Γi_q = Γ_q[iΓ+i] - Γ_q[iΓ+i-1]
                Γi_r = Γ_r[iΓ+i] - Γ_r[iΓ+i-1]
            end

            # bound vortex vector
            Δs = top_vector(receiving[I])

            tmp = cross(Vi, Δs)

            Fbi = ref.rho*Γi*tmp

            Fbi_a = ref.rho*(Γi_a*tmp + Γi*cross(Vi_a, Δs))
            Fbi_b = ref.rho*(Γi_b*tmp + Γi*cross(Vi_b, Δs))
            Fbi_p = ref.rho*(Γi_p*tmp + Γi*cross(Vi_p, Δs))
            Fbi_q = ref.rho*(Γi_q*tmp + Γi*cross(Vi_q, Δs))
            Fbi_r = ref.rho*(Γi_r*tmp + Γi*cross(Vi_r, Δs))

            if !isnothing(dΓdt)
                # unsteady part of Kutta-Joukowski theorem

                #TODO: decide whether to divide by perpindicular velocity like
                # Drela does in ASWING?

                dΓdti = I[1] == 1 ? dΓdt[iΓ+i] : (dΓdt[iΓ+i] + dΓdt[iΓ+i-1])/2
                c = receiving[I].chord
                Fbi += ref.rho*dΓdti*c*tmp

            end

            # --- Calculate forces for the left bound vortex --- #

            rc = left_center(receiving[I])

            # freestream velocity
            Vfs, dVfs = freestream_velocity_derivatives(fs)
            Veff = Vfs
            Veff_a, Veff_b = dVfs

            # rotational velocity
            Vrot, dVrot = rotational_velocity_derivatives(rc, fs, ref)
            Veff += Vrot
            Veff_p, Veff_q, Veff_r = dVrot

            # additional velocity field
            if !isnothing(additional_velocity)
                Veff += additional_velocity(rc)
            end

            # velocity due to surface motion
            if !isnothing(Vv)
                Veff += Vv[isurf][I[1], I[2]]
            end

            # NOTE: We don't include induced velocity in the effective velocity
            # for the vertical segments because its influence is likely negligible
            # once we take the cross product with the bound vortex vector. This
            # is also assumed in AVL. This could change in the future.

            # steady part of Kutta-Joukowski theorem
            Γli = Γ[iΓ+i]

            Γli_a = Γ_a[iΓ+i]
            Γli_b = Γ_b[iΓ+i]
            Γli_p = Γ_p[iΓ+i]
            Γli_q = Γ_q[iΓ+i]
            Γli_r = Γ_r[iΓ+i]

            Δs = left_vector(receiving[I])

            tmp = cross(Veff, Δs)

            Fbli = ref.rho*Γli*tmp

            Fbli_a = ref.rho*(Γli_a*tmp + Γli*cross(Veff_a, Δs))
            Fbli_b = ref.rho*(Γli_b*tmp + Γli*cross(Veff_b, Δs))
            Fbli_p = ref.rho*(Γli_p*tmp + Γli*cross(Veff_p, Δs))
            Fbli_q = ref.rho*(Γli_q*tmp + Γli*cross(Veff_q, Δs))
            Fbli_r = ref.rho*(Γli_r*tmp + Γli*cross(Veff_r, Δs))

            # --- Calculate forces on the right bound vortex --- #

            rc = right_center(receiving[I])

            # freestream velocity
            Vfs, dVfs = freestream_velocity_derivatives(fs)
            Veff = Vfs
            Veff_a, Veff_b = dVfs

            # rotational velocity
            Vrot, dVrot = rotational_velocity_derivatives(rc, fs, ref)
            Veff += Vrot
            Veff_p, Veff_q, Veff_r = dVrot

            # additional velocity field
            if !isnothing(additional_velocity)
                Veff += additional_velocity(rc)
            end

            # velocity due to surface motion
            if !isnothing(Vv)
                Veff += Vv[isurf][I[1], I[2]+1]
            end

            # NOTE: We don't include induced velocity in the effective velocity
            # for the vertical segments because its influence is likely negligible
            # once we take the cross product with the bound vortex vector. This
            # is also assumed in AVL. This could change in the future.

            # steady part of Kutta-Joukowski theorem
            Γri = Γ[iΓ+i]

            Γri_a = Γ_a[iΓ+i]
            Γri_b = Γ_b[iΓ+i]
            Γri_p = Γ_p[iΓ+i]
            Γri_q = Γ_q[iΓ+i]
            Γri_r = Γ_r[iΓ+i]

            Δs = right_vector(receiving[I])

            tmp = cross(Veff, Δs)

            Fbri = ref.rho*Γri*tmp

            Fbri_a = ref.rho*(Γri_a*tmp + Γri*cross(Veff_a, Δs))
            Fbri_b = ref.rho*(Γri_b*tmp + Γri*cross(Veff_b, Δs))
            Fbri_p = ref.rho*(Γri_p*tmp + Γri*cross(Veff_p, Δs))
            Fbri_q = ref.rho*(Γri_q*tmp + Γri*cross(Veff_q, Δs))
            Fbri_r = ref.rho*(Γri_r*tmp + Γri*cross(Veff_r, Δs))

            # store panel circulation, velocity, and forces
            q = 1/2*ref.rho*ref.V^2

            props[isurf][I] = PanelProperties(Γ[iΓ+i]/ref.V, Vi/ref.V, Fbi,
                Fbli, Fbri, V_streamwise)

            props_a[isurf][I] = PanelProperties(Γ_a[iΓ+i]/ref.V, Vi_a/ref.V, Fbi_a/(q*ref.S),
                Fbli_a/(q*ref.S), Fbri_a/(q*ref.S), V_streamwise)
            props_b[isurf][I] = PanelProperties(Γ_b[iΓ+i]/ref.V, Vi_b/ref.V, Fbi_b/(q*ref.S),
                Fbli_b/(q*ref.S), Fbri_b/(q*ref.S), V_streamwise)
            props_p[isurf][I] = PanelProperties(Γ_p[iΓ+i]/ref.V, Vi_p/ref.V, Fbi_p/(q*ref.S),
                Fbli_p/(q*ref.S), Fbri_p/(q*ref.S), V_streamwise)
            props_q[isurf][I] = PanelProperties(Γ_q[iΓ+i]/ref.V, Vi_q/ref.V, Fbi_q/(q*ref.S),
                Fbli_q/(q*ref.S), Fbri_q/(q*ref.S), V_streamwise)
            props_r[isurf][I] = PanelProperties(Γ_r[iΓ+i]/ref.V, Vi_r/ref.V, Fbi_r/(q*ref.S),
                Fbli_r/(q*ref.S), Fbri_r/(q*ref.S), V_streamwise)
        end

        # increment Γ index for receiving panels
        iΓ += nr
    end

    return props, dprops
end

"""
    body_forces(system; kwargs...)

Return the body force coefficients given the panel properties for `surfaces`

Note that this function assumes that a near-field analysis has already been
performed to obtain the panel forces.

# Arguments
 - `system`: Object of type [`System`](@ref) which holds system properties

# Keyword Arguments
 - `frame`: frame in which to return `CF` and `CM`, options are [`Body()`](@ref) (default),
   [`Stability()`](@ref), and [`Wind()`](@ref)`
"""
function body_forces(system::System{TF}; frame = Body()) where TF

    @assert system.near_field_analysis[] "Near field analysis required"

    # unpack parameters stored in `system`
    surfaces = system.surfaces # surface panels defining each surface
    properties = system.properties # panel properties
    ref = system.reference[] # reference parameters
    fs = system.freestream[] # freestream parameters
    symmetric = system.symmetric # symmetric flag for each surface

    return body_forces(surfaces, properties, ref, fs, symmetric, frame)
end

"""
    body_forces(surfaces, properties, reference, freestream, symmetric; kwargs...)

Return the body force coefficients given the panel properties for `surfaces`

Note that this function assumes that a near-field analysis has already been
performed to obtain the panel forces.

# Arguments:
 - `surfaces`: Collection of surfaces, where each surface is represented by a
    matrix of surface panels (see [`SurfacePanel`](@ref)) of shape (nc, ns)
    where `nc` is the number of chordwise panels and `ns` is the number of
    spanwise panels
 - `properties`: Surface properties for each surface, where surface
    properties for each surface are represented by a matrix of panel properties
    (see [`PanelProperties`](@ref)) of shape (nc, ns) where `nc` is the number
    of chordwise panels and `ns` is the number of spanwise panels
 - `reference`: Reference parameters (see [`Reference`](@ref))
 - `freestream`: Freestream parameters (see [`Freestream`]@ref)
 - `symmetric`: (required) Flag for each surface indicating whether a mirror image
   (across the X-Z plane) was used when calculating induced velocities
 - `frame`: frame in which to return `CF` and `CM`, options are [`Body()`](@ref) (default),
   [`Stability()`](@ref), and [`Wind()`](@ref)
"""
function body_forces(surfaces, properties, ref, fs, symmetric, frame)

    TF = eltype(eltype(eltype(properties)))

    # initialize body force coefficients
    CF = @SVector zeros(TF, 3)
    CM = @SVector zeros(TF, 3)

    # loop through all surfaces
    for isurf = 1:length(surfaces)

        # initialize surface contribution to body force coefficients
        CFi = @SVector zeros(TF, 3)
        CMi = @SVector zeros(TF, 3)

        # loop through all panels on this surface
        for i = 1:length(surfaces[isurf])

            # top bound vortex
            rc = top_center(surfaces[isurf][i])
            Δr = rc - ref.r
            cf = properties[isurf][i].cfb
            CFi += cf
            CMi += cross(Δr, cf)

            # left bound vortex
            rc = left_center(surfaces[isurf][i])
            Δr = rc - ref.r
            cf = properties[isurf][i].cfl
            CFi += cf
            CMi += cross(Δr, cf)

            # right bound vortex
            rc = right_center(surfaces[isurf][i])
            Δr = rc - ref.r
            cf = properties[isurf][i].cfr
            CFi += cf
            CMi += cross(Δr, cf)
        end

        # adjust forces from this surface to account for symmetry
        if symmetric[isurf]
            CFi = SVector(2*CFi[1], 0.0, 2*CFi[3])
            CMi = SVector(0.0, 2*CMi[2], 0.0)
        end

        # add to body forces
        CF += CFi
        CM += CMi

    end

    # add reference length in moment normalization
    reference_length = SVector(ref.b, ref.c, ref.b)
    CM = CM ./ reference_length

    # positive Mx corresponds to negative roll, and positive Mz corresponds to negative yaw
    convention_change = [-1.0, 1.0, -1.0]
    CM = CM .* convention_change

    # switch to specified frame
    CF, CM = body_to_frame(CF, CM, ref, fs, frame)

    return CF, CM
end

"""
    body_forces_derivatives(system)

Return the body force coefficients for the `system` and their derivatives with
respect to the freestream variables

Note that this function assumes that a near-field analysis has already been
performed to obtain the panel forces.

# Arguments:
 - `system`: Object of type `System` which holds system properties
"""
function body_forces_derivatives(system::System)

    # float number type
    TF = eltype(system)

    @assert system.near_field_analysis[] "Near field analysis required"
    @assert system.derivatives[] "Derivative computations required"

    # unpack parameters stored in `system`
    surfaces = system.surfaces # surface panels defining each surface
    ref = system.reference[] # reference parameters
    fs = system.freestream[] # freestream parameters
    symmetric = system.symmetric # symmetric flag for each surface
    properties = system.properties
    props_a, props_b, props_p, props_q, props_r = system.dproperties

    # initialize body force coefficients
    CF = @SVector zeros(TF, 3)
    CM = @SVector zeros(TF, 3)

    CF_a = @SVector zeros(TF, 3)
    CF_b = @SVector zeros(TF, 3)
    CF_p = @SVector zeros(TF, 3)
    CF_q = @SVector zeros(TF, 3)
    CF_r = @SVector zeros(TF, 3)

    CM_a = @SVector zeros(TF, 3)
    CM_b = @SVector zeros(TF, 3)
    CM_p = @SVector zeros(TF, 3)
    CM_q = @SVector zeros(TF, 3)
    CM_r = @SVector zeros(TF, 3)

    # loop through all surfaces
    for isurf = 1:length(surfaces)

        # initialize surface contribution to body force coefficients
        CFi = @SVector zeros(TF, 3)
        CMi = @SVector zeros(TF, 3)

        CFi_a = @SVector zeros(TF, 3)
        CFi_b = @SVector zeros(TF, 3)
        CFi_p = @SVector zeros(TF, 3)
        CFi_q = @SVector zeros(TF, 3)
        CFi_r = @SVector zeros(TF, 3)

        CMi_a = @SVector zeros(TF, 3)
        CMi_b = @SVector zeros(TF, 3)
        CMi_p = @SVector zeros(TF, 3)
        CMi_q = @SVector zeros(TF, 3)
        CMi_r = @SVector zeros(TF, 3)

        for i = 1:length(surfaces[isurf])

            # top bound vortex
            rc = top_center(surfaces[isurf][i])
            Δr = rc - ref.r
            cf = properties[isurf][i].cfb
            CFi += cf
            CMi += cross(Δr, cf)

            cf_a = props_a[isurf][i].cfb
            CFi_a += cf_a
            CMi_a += cross(Δr, cf_a)

            cf_b = props_b[isurf][i].cfb
            CFi_b += cf_b
            CMi_b += cross(Δr, cf_b)

            cf_p = props_p[isurf][i].cfb
            CFi_p += cf_p
            CMi_p += cross(Δr, cf_p)

            cf_q = props_q[isurf][i].cfb
            CFi_q += cf_q
            CMi_q += cross(Δr, cf_q)

            cf_r = props_r[isurf][i].cfb
            CFi_r += cf_r
            CMi_r += cross(Δr, cf_r)

            # left bound vortex
            rc = left_center(surfaces[isurf][i])
            Δr = rc - ref.r
            cf = properties[isurf][i].cfl
            CFi += cf
            CMi += cross(Δr, cf)

            cf_a = props_a[isurf][i].cfl
            CFi_a += cf_a
            CMi_a += cross(Δr, cf_a)

            cf_b = props_b[isurf][i].cfl
            CFi_b += cf_b
            CMi_b += cross(Δr, cf_b)

            cf_p = props_p[isurf][i].cfl
            CFi_p += cf_p
            CMi_p += cross(Δr, cf_p)

            cf_q = props_q[isurf][i].cfl
            CFi_q += cf_q
            CMi_q += cross(Δr, cf_q)

            cf_r = props_r[isurf][i].cfl
            CFi_r += cf_r
            CMi_r += cross(Δr, cf_r)

            # right bound vortex
            rc = right_center(surfaces[isurf][i])
            Δr = rc - ref.r
            cf = properties[isurf][i].cfr
            CFi += cf
            CMi += cross(Δr, cf)

            cf_a = props_a[isurf][i].cfr
            CFi_a += cf_a
            CMi_a += cross(Δr, cf_a)

            cf_b = props_b[isurf][i].cfr
            CFi_b += cf_b
            CMi_b += cross(Δr, cf_b)

            cf_p = props_p[isurf][i].cfr
            CFi_p += cf_p
            CMi_p += cross(Δr, cf_p)

            cf_q = props_q[isurf][i].cfr
            CFi_q += cf_q
            CMi_q += cross(Δr, cf_q)

            cf_r = props_r[isurf][i].cfr
            CFi_r += cf_r
            CMi_r += cross(Δr, cf_r)
        end

        # adjust forces from this surface to account for symmetry
        if symmetric[isurf]
            CFi = SVector(2*CFi[1], 0.0, 2*CFi[3])
            CMi = SVector(0.0, 2*CMi[2], 0.0)

            CFi_a = SVector(2*CFi_a[1], 0.0, 2*CFi_a[3])
            CFi_b = SVector(2*CFi_b[1], 0.0, 2*CFi_b[3])
            CFi_p = SVector(2*CFi_p[1], 0.0, 2*CFi_p[3])
            CFi_q = SVector(2*CFi_q[1], 0.0, 2*CFi_q[3])
            CFi_r = SVector(2*CFi_r[1], 0.0, 2*CFi_r[3])

            CMi_a = SVector(0.0, 2*CMi_a[2], 0.0)
            CMi_b = SVector(0.0, 2*CMi_b[2], 0.0)
            CMi_p = SVector(0.0, 2*CMi_p[2], 0.0)
            CMi_q = SVector(0.0, 2*CMi_q[2], 0.0)
            CMi_r = SVector(0.0, 2*CMi_r[2], 0.0)
        end

        # add to body forces
        CF += CFi
        CM += CMi

        CF_a += CFi_a
        CF_b += CFi_b
        CF_p += CFi_p
        CF_q += CFi_q
        CF_r += CFi_r

        CM_a += CMi_a
        CM_b += CMi_b
        CM_p += CMi_p
        CM_q += CMi_q
        CM_r += CMi_r

    end

    # add reference length in moment normalization
    reference_length = SVector(ref.b, ref.c, ref.b)
    CM = CM ./ reference_length

    CM_a = CM_a ./ reference_length
    CM_b = CM_b ./ reference_length
    CM_p = CM_p ./ reference_length
    CM_q = CM_q ./ reference_length
    CM_r = CM_r ./ reference_length

    # positive Mx corresponds to negative roll, and positive Mz corresponds to negative yaw
    convention_change = [-1.0, 1.0, -1.0]
    CM_a = CM_a .* convention_change
    CM_b = CM_b .* convention_change
    CM_p = CM_p .* convention_change
    CM_q = CM_q .* convention_change
    CM_r = CM_r .* convention_change

    # pack up derivatives
    dCF = (CF_a, CF_b, CF_p, CF_q, CF_r)
    dCM = (CM_a, CM_b, CM_p, CM_q, CM_r)

    return CF, CM, dCF, dCM
end

"""
    body_forces_history(system, surface_history, property_history; frame=Body())

Return the body force coefficients `CF`, `CM` at each time step in `property_history`.

# Arguments:
 - `system`: Object of type [`System`](@ref) which holds system properties
 - `surface_history`: Vector of surfaces at each time step, where each surface is
    represented by a matrix of surface panels (see [`SurfacePanel`](@ref)) of shape
    (nc, ns) where `nc` is the number of chordwise panels and `ns` is the number
    of spanwise panels
 - `property_history`: Vector of surface properties for each surface at each
    time step, where surface properties are represented by a matrix of panel
    properties (see [`PanelProperties`](@ref)) of shape (nc, ns) where `nc` is
    the number of chordwise panels and `ns` is the number of spanwise panels

# Keyword Arguments
 - `frame`: frame in which to return `CF` and `CM`, options are [`Body()`](@ref) (default),
   [`Stability()`](@ref), and [`Wind()`](@ref)`
"""
function body_forces_history(system, surface_history::AbstractVector{<:AbstractVector{<:AbstractMatrix}},
    property_history; frame=Body())

    # unpack system parameters
    symmetric = system.symmetric
    ref = system.reference[]
    fs = system.freestream[]

    # float type
    TF = eltype(system)

    # number of time steps
    nt = length(property_history)

    # convert single freestream input to vector
    if isa(fs, Freestream)
        fs = fill(fs, nt)
    end

    # initialize time history coefficients
    CF = Vector{SVector{3, TF}}(undef, nt)
    CM = Vector{SVector{3, TF}}(undef, nt)

    # populate time history coefficients
    for it = 1:nt
        CF[it], CM[it] = body_forces(surface_history[it], property_history[it],
            ref, fs[it], symmetric, frame)
    end

    return CF, CM
end

"""
    lifting_line_coefficients(system, r, c; frame=Body())

Return the force and moment coefficients (per unit span) for each spanwise segment
of a lifting line representation of the geometry.

This function requires that a near-field analysis has been performed on `system`
to obtain panel forces.

# Arguments
 - `system`: Object of type [`System`](@ref) that holds precalculated
    system properties.
 - `r`: Vector with length equal to the number of surfaces, with each element
    being a matrix with size (3, ns+1) which contains the x, y, and z coordinates
    of the resulting lifting line coordinates
 - `c`: Vector with length equal to the number of surfaces, with each element
    being a vector of length `ns+1` which contains the chord lengths at each
    lifting line coordinate.

# Keyword Arguments
 - `frame`: frame in which to return `cf` and `cm`, possible options are
    [`Body()`](@ref) (default), [`Stability()`](@ref), and [`Wind()`](@ref)`

# Return Arguments:
 - `cf`: Vector with length equal to the number of surfaces, with each element
    being a matrix with size (3, ns) which contains the x, y, and z direction
    force coefficients (per unit span) for each spanwise segment.
 - `cm`: Vector with length equal to the number of surfaces, with each element
    being a matrix with size (3, ns) which contains the x, y, and z direction
    moment coefficients (per unit span) for each spanwise segment.
"""
function lifting_line_coefficients(system, r, c; frame=Body())
    TF = promote_type(eltype(system), eltype(eltype(r)), eltype(eltype(c)))
    nsurf = length(system.surfaces)
    cf = Vector{Matrix{TF}}(undef, nsurf)
    cm = Vector{Matrix{TF}}(undef, nsurf)
    for isurf = 1:nsurf
        ns = size(system.surfaces[isurf], 2)
        cf[isurf] = Matrix{TF}(undef, 3, ns)
        cm[isurf] = Matrix{TF}(undef, 3, ns)
    end
    return lifting_line_coefficients!(cf, cm, system, r, c; frame)
end

function lifting_line_coefficients(system; frame=Body(), xc = 0.25)
    r, c = lifting_line_geometry(system.grids, xc)
    return lifting_line_coefficients(system, r, c; frame)
end

"""
    lifting_line_coefficients!(cf, cm, system, r, c; frame=Body())

In-place version of [`lifting_line_coefficients`](@ref)
"""
function lifting_line_coefficients!(cf, cm, system, r, c; frame=Body())

    # number of surfaces
    nsurf = length(system.surfaces)

    # check that a near field analysis has been performed
    @assert system.near_field_analysis[] "Near field analysis required"

    # extract reference properties
    ref = system.reference[]
    fs = system.freestream[]

    # iterate through each lifting surface
    for isurf = 1:nsurf
        nc, ns = size(system.surfaces[isurf])
        # extract current surface panels and panel properties
        panels = system.surfaces[isurf]
        properties = system.properties[isurf]
        # loop through each chordwise set of panels
        for j = 1:ns
            # calculate segment length
            rls = SVector(r[isurf][1,j], r[isurf][2,j], r[isurf][3,j])
            rrs = SVector(r[isurf][1,j+1], r[isurf][2,j+1], r[isurf][3,j+1])
            ds = norm(rrs - rls)
            # calculate reference location
            rs = (rls + rrs)/2
            # calculate reference chord
            cs = (c[isurf][j] + c[isurf][j+1])/2
            # calculate section force and moment coefficients
            cfj = @SVector zeros(eltype(cf[isurf]), 3)
            cmj = @SVector zeros(eltype(cm[isurf]), 3)
            for i = 1:nc
                # add influence of bound vortex
                rb = top_center(panels[i,j])
                #rb = controlpoint(panels[i,j])
                #controlpoint
                cfb = properties[i,j].cfb
                cfj += cfb
                cmj += cross(rb - rs, cfb)
                # add influence of left vortex leg
                rl = left_center(panels[i,j])
                #rl = controlpoint(panels[i,j])
                cfl = properties[i,j].cfl
                cfj += cfl
                cmj += cross(rl - rs, cfl)
                # add influence of right vortex leg
                rr = right_center(panels[i,j])
                #rr = controlpoint(panels[i,j])
                cfr = properties[i,j].cfr
                cfj += cfr
                cmj += cross(rr - rs, cfr)
            end
            # update normalization
            cfj *= ref.S/(ds*cs)
            cmj *= ref.S/(ds*cs^2)
            # change coordinate frame
            cfj, cmj = body_to_frame(cfj, cmj, ref, fs, frame)
            # save coefficients
            cf[isurf][:,j] = cfj
            cm[isurf][:,j] = cmj
        end
    end

    return cf, cm
end

"""
    body_to_frame(CF, CM, reference, freestream, frame)

Transform the coefficients `CF` and `CM` from the body frame to the frame
specified in `frame`
"""
body_to_frame

body_to_frame(CF, CM, ref, fs, ::Body) = CF, CM

function body_to_frame(CF, CM, ref, fs, ::Stability)
    R = body_to_stability(fs)
    return R*CF, R*CM
end

function body_to_frame(CF, CM, ref, fs, ::Wind)
    # remove reference lengths
    reflen = SVector(ref.b, ref.c, ref.b)
    CM = CM .* reflen

    # rotate
    R = body_to_wind(fs)
    CF = R*CF
    CM = R*CM

    # reapply reference lengths
    CM = CM ./ reflen

    return CF, CM
end
