# The bridge panel's rear edge is opposite the first wake panel's front edge.
# Use the wake kernel for both terms, retaining their respective strengths.
function match_wake_interface_core!(AIC,surfaces,wakes,nwake,surface_id,wake_finite_core,
    symmetric,trailing_vortices,interaction_id,interaction)

    offsets = cumsum(vcat(0,length.(surfaces)))
    for isend in eachindex(surfaces)
        (nwake[isend] > 0 && !trailing_vortices[isend]) || continue
        sending = surfaces[isend]
        nc,ns = size(sending)
        for irecv in eachindex(surfaces)
            (interaction || interaction_id[irecv] == interaction_id[isend]) || continue
            sameSurface = surface_id[irecv] == surface_id[isend]
            boundCore = !sameSurface
            wakeCore = wake_finite_core[isend] || !sameSurface
            receiving = surfaces[irecv]
            for j = 1:ns
                firstWake = wakes[isend][1,j]
                oldRadius = get_core_size(sending[end,j])
                newRadius = get_core_size(firstWake)
                boundCore == wakeCore && (!wakeCore || oldRadius == newRadius) && continue
                left,right = top_left(firstWake),top_right(firstWake)
                column = offsets[isend] + nc * j
                for i in eachindex(receiving)
                    cp = controlpoint(receiving[i])
                    r1,r2 = cp - right,cp - left
                    difference = bound_induced_velocity(r1,r2,wakeCore,newRadius) -
                                 bound_induced_velocity(r1,r2,boundCore,oldRadius)
                    if symmetric[isend]
                        r1,r2 = cp - flipy(left),cp - flipy(right)
                        difference += bound_induced_velocity(r1,r2,wakeCore,newRadius) -
                                      bound_induced_velocity(r1,r2,boundCore,oldRadius)
                    end
                    AIC[offsets[irecv]+i,column] += dot(difference,normal(receiving[i]))
                end
            end
        end
    end
    return AIC
end
