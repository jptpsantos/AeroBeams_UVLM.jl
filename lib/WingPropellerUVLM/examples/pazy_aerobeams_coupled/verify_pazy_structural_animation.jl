# Render a short preview from a saved result, without rerunning the simulation.
# Usage: julia --project=. verify_pazy_structural_animation.jl path/to/result.jls
using Serialization, Test
include(joinpath(@__DIR__, "PazyWingUVLMCoupling.jl"))
include(joinpath(@__DIR__, "PazyWingUVLMVisualization.jl"))

function verify_animation(result)
    # Work on a copy: preserve the stored solution and its full frame history.
    preview = deepcopy(result)
    problem = preview.structural_problem
    bcs = problem.model.BCs
    surface = problem.model.beams[1].aeroSurface
    element_aero = [e.aero for e in problem.model.elements]
    bc_types = [copy(bc.types) for bc in bcs]
    bc_values = [copy(bc.values) for bc in bcs]
    render_env = (get(ENV,"MKL_THREADING_LAYER",nothing),get(ENV,"MKL_NUM_THREADS",nothing))
    # The investigation checkpoint kept only initial structural history. Append
    # its actual final solved state for this preview (no interpolation).
    if last(problem.savedTimeVector) < last(preview.time)
        AeroBeams.save_time_step_data!(problem, last(preview.time))
    end
    preview = merge(preview, (; animation_stride=max(1,length(problem.savedTimeVector)-1)))
    folder = joinpath(@__DIR__, "output", "structural_style_preview")
    path = save_structural_animation(preview;
        output_path=joinpath(folder,"pazy_structure_preview.gif"), fps=1,
        camera=(45,45), surface_alpha=0.5, elastic_axis_linewidth=2.5,
        show_clamped_root=true, show_force_vectors=false, show_moment_vectors=false)
    Plots.savefig(joinpath(folder,"pazy_structure_preview.png"))
    @testset "Structural animation preserves the solution" begin
        @test isfile(path)
        @test problem.model.BCs === bcs
        @test problem.model.beams[1].aeroSurface === surface
        @test all(e.aero === original for (e,original) in zip(problem.model.elements,element_aero))
        @test [bc.types for bc in bcs] == bc_types
        @test [bc.values for bc in bcs] == bc_values
        @test problem.x == result.structural_problem.x
        @test problem.savedTimeVector[end] == last(result.time)
        @test (get(ENV,"MKL_THREADING_LAYER",nothing),get(ENV,"MKL_NUM_THREADS",nothing)) == render_env
    end
    println("Preview only (first/last stored states), not a new simulation: ",folder)
end

data = deserialize(only(ARGS))
Base.invokelatest(verify_animation,data)
