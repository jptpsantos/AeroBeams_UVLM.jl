# Short end-to-end smoke test for all independent integrator/coupling choices.
# Run with --project=lib/WingPropellerUVLM from the repository root.
using Test
using LinearAlgebra
using WingPropellerUVLM

include(joinpath(@__DIR__, "..", "src", "ChangAeroelastic.jl"))
using .ChangAeroelastic

BLAS.set_num_threads(1)

function combination_smoke_defaults(
    time_integrator,
    coupling_scheme,
    label;
    azimuth_step_deg = 90.0,
)
    defaults = chang_case_defaults()
    return merge(defaults, (;
        wing = merge(defaults.wing, (;
            spanwise_panels = 2,
            chordwise_panels = 1,
        )),
        propeller = merge(defaults.propeller, (;
            blades = 2,
            radial_panels = 1,
            chordwise_panels = 1,
        )),
        simulation = merge(defaults.simulation, (;
            azimuth_step_deg,
            end_time_s = 0.05,
            interaction_on = false,
        )),
        wake = merge(defaults.wake, (;
            maximum_rows_wing = 16,
            maximum_rows_propeller = 16,
        )),
        excitation = merge(defaults.excitation, (;
            trim_revolutions = 0.25,
            trim_average_revolutions = 0.25,
            impulse_duration_s = 0.02,
            impulse_magnitude_nm = 10.0,
        )),
        integration = merge(defaults.integration, (;
            time_integrator,
            state_norm_limit = 1.0e6,
        )),
        coupling = merge(defaults.coupling, (;
            scheme = coupling_scheme,
            maximum_iterations = 30,
            state_tolerance = 1.0e-5,
            load_tolerance = 1.0e-3,
            equilibrium_tolerance = 1.0e-8,
            coupled_equilibrium_tolerance = 1.0e-3,
            relaxation = 1.0,
            verbose = false,
        )),
        output = merge(defaults.output, (;
            directory = joinpath(@__DIR__, "..", "output", "solver_combinations"),
            label,
            plot_results = false,
            animate_wake = false,
        )),
    ))
end

@testset "Chang solver combinations" begin
    runs = Dict{Tuple{Symbol,Symbol},Any}()
    for time_integrator in STRUCTURAL_TIME_INTEGRATORS
        for coupling_scheme in AEROELASTIC_COUPLING_SCHEMES
            label = "smoke_$(time_integrator)_$(coupling_scheme)"
            defaults = combination_smoke_defaults(
                time_integrator, coupling_scheme, label,
            )
            config = load_chang_configuration(defaults; env = Dict{String,String}())
            run = run_chang(config)
            runs[(time_integrator, coupling_scheme)] = run
            solution = run.solution
            active_steps = 1:solution.last_step

            @test solution.time_integrator == time_integrator
            @test solution.coupling_scheme == coupling_scheme
            @test solution.last_step == length(run.model.parameters.dt)
            @test all(solution.coupling_converged[active_steps])
            @test all(
                all(isfinite, state)
                for state in solution.displacement_history[1:(solution.last_step + 1)]
            )
            @test all(run.workspace.iwake .== solution.last_step)
            if coupling_scheme == :loose_explicit
                @test all(solution.coupling_iterations[active_steps] .== 1)
            else
                @test any(solution.coupling_iterations[active_steps] .> 1)
                @test maximum(solution.coupling_state_residual[active_steps]) <=
                    config.coupling.state_tolerance
            end
        end
    end

    # Compare the loose/implicit histories again at half the azimuth/time
    # increment. This is diagnostic rather than a monotonicity assertion: the
    # free-wake temporal discretization and number of trim-average samples also
    # change with azimuth step, so a five-step smoke transient is not a formal
    # time-convergence or damping study.
    refined_runs = Dict{Symbol,Any}()
    for coupling_scheme in AEROELASTIC_COUPLING_SCHEMES
        label = "refined_newmark_beta_$(coupling_scheme)"
        defaults = combination_smoke_defaults(
            :newmark_beta, coupling_scheme, label; azimuth_step_deg = 45.0,
        )
        config = load_chang_configuration(defaults; env = Dict{String,String}())
        refined_runs[coupling_scheme] = run_chang(config)
    end
    coarse_explicit = runs[(:newmark_beta, :loose_explicit)].solution
    coarse_implicit = runs[(:newmark_beta, :implicit_predictor_corrector)].solution
    refined_explicit = refined_runs[:loose_explicit].solution
    refined_implicit = refined_runs[:implicit_predictor_corrector].solution
    coarse_difference = maximum(
        norm(coarse_explicit.displacement_history[index] .-
            coarse_implicit.displacement_history[index])
        for index in 1:(coarse_explicit.last_step + 1)
    )
    refined_difference = maximum(
        norm(refined_explicit.displacement_history[2 * index - 1] .-
            refined_implicit.displacement_history[2 * index - 1])
        for index in 1:(coarse_explicit.last_step + 1)
    )
    @test isfinite(coarse_difference) && isfinite(refined_difference)
    println(
        "Newmark loose/implicit full-state difference: coarse=$coarse_difference, " *
        "half-step=$refined_difference",
    )
end
