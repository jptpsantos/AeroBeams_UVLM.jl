# Editable case. This file constructs settings only; the runner performs work.
defaults = ChangAeroBeams.chang_aerobeams_defaults()

# Keep mean-equilibrium initialization for physical loaded-wing analyses.
# For short numerical checks choose mode=:unloaded and use a small duration.
case = merge(
    defaults,
    (
        reference=merge(
            defaults.reference,
            (
                propeller=merge(defaults.reference.propeller,(attachment_eta=[0.65],)),
                simulation=merge(
                    defaults.reference.simulation,
                    (
                        freestream_speed_mps=65.0,end_time_s=0.25,interaction_on=true)
                )
            )
        ),
        solver=merge(defaults.solver,(support=:reduced,)),
        initialization=merge(defaults.initialization,(mode=:mean_equilibrium,)),
        pylon=merge(defaults.pylon,(elements=4,)),
        output=merge(defaults.output,(plot_results=true,))
    )
)
case
