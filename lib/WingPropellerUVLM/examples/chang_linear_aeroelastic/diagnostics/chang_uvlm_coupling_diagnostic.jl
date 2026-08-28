# Diagnostic wrapper around the production Chang coupling adapter.
# It can replace only the near-field force reconstruction while leaving the
# circulation, wake, geometry, structure, and time integrator unchanged.

import WingPropellerUVLM

const PRODUCTION_COUPLING_FILE = normpath(joinpath(@__DIR__, "..", "chang_uvlm_coupling.jl"))
const TEST_FORCE_MODEL = Symbol(lowercase(strip(get(ENV, "CHANG_TEST_FORCE_MODEL", "imperial"))))
TEST_FORCE_MODEL in (:imperial, :legacy) ||
    error("CHANG_TEST_FORCE_MODEL must be imperial or legacy")

function diagnostic_select_nearfield_forces!(system)
    TEST_FORCE_MODEL == :imperial && return nothing

    legacy_properties = deepcopy(system.properties)
    _, chord_forces, span_forces, unsteady_forces =
        WingPropellerUVLM.legacy_near_field_forces!(
            legacy_properties,
            system.surfaces,
            system.wakes,
            system.reference[],
            system.freestream[],
            system.Γ;
            dΓdt = system.dΓdt,
            additional_velocity = nothing,
            Vh = system.Vh,
            Vv = system.Vv,
            symmetric = system.symmetric,
            nwake = system.nwake,
            surface_id = system.surface_id,
            wake_finite_core = system.wake_finite_core,
            wake_shedding_locations = system.wake_shedding_locations,
            trailing_vortices = system.trailing_vortices,
            xhat = system.xhat[],
            interaction_id = surface_interaction_id,
            interaction = INTERACTION_ON,
        )

    for surface_index in eachindex(system.surfaces)
        system.chord_seg_forces[surface_index] .= chord_forces[surface_index]
        system.span_seg_forces[surface_index] .= span_forces[surface_index]
        system.unsteady_forces[surface_index] .= unsteady_forces[surface_index]
    end
    return nothing
end

coupling_source = read(PRODUCTION_COUPLING_FILE, String)
coupling_source = replace(
    coupling_source,
    "    # Convert the aerodynamic trial into a load vector for the structural\n" *
    "    # corrector. The caller decides whether to subtract the trim baseline.\n" *
    "    return assemble_structural_aero_load!(" =>
    "    # Diagnostic-only selection of the force reconstruction.\n" *
    "    diagnostic_select_nearfield_forces!(system)\n" *
    "    # Convert the aerodynamic trial into a load vector for the structural\n" *
    "    # corrector. The caller decides whether to subtract the trim baseline.\n" *
    "    return assemble_structural_aero_load!(",
)

Base.include_string(@__MODULE__, coupling_source, PRODUCTION_COUPLING_FILE)
