using Test
using Serialization
include(joinpath(@__DIR__,"..","PazyWingUVLMStatic.jl"))

@testset "Native Pazy options and swept UVLM geometry" begin
    cases = ((;),(; useAePW3Matrices=true),(; mesh=:refined),
        (; withWingtipRod=false),(; withWingtipRod=false,mesh=:refined),
        (; rootTorsionCorr=false),(; beamWingtipRod=true),
        (; sweepAngle=deg2rad(10)),(; sweepAngle=deg2rad(20),mesh=:refined),
        (; sweepAngle=deg2rad(10),sweptTipMassConfig=:LE,sweptTipMassInertiaCorrection=true),
        (; sweepAngle=deg2rad(20),sweptTipMassConfig=:TE),
        (; sweepAngle=deg2rad(10),beamWingtipRod=true,sweptTipMassConfig=:LE),
        (; sweepAngle=deg2rad(10),beamWingtipRod=true,sweptTipMassConfig=:TE))
    for options in cases
        @testset "$options" begin
            pazy = create_pazy_uvlm_structure(; options...)
            native,ne,L,c,spar = AeroBeams.create_Pazy(; pazy.parameters...,
                upright=true,pitchAngle=0.0,gravityVector=zeros(3),hasTipCorrection=false)
            @test (pazy.n_elements,pazy.span,pazy.chord,pazy.sparFraction) == (ne,L,c,spar)
            @test length(pazy.model.beams) == length(native.beams)
            @test length(pazy.model.connections) == length(native.connections)
            for (beam,referenceBeam) in zip(pazy.model.beams,native.beams)
                @test beam.sectional.S == referenceBeam.sectional.S
                @test beam.sectional.I == referenceBeam.sectional.I
                @test length(beam.attachments.springs) == length(referenceBeam.attachments.springs)
                @test length(beam.attachments.pointInertias) == length(referenceBeam.attachments.pointInertias)
                @test isnothing(beam.attachments.surface)
                @test all(isnothing(e.aero.surface) for e in beam.elements)
            end
            @test length(pazy.model.BCs) == ne + 2
            inertia = pazy.model.assemblyData.initialInertiaProperties
            referenceInertia = native.assemblyData.initialInertiaProperties
            @test inertia.mass ≈ referenceInertia.mass
            @test inertia.centerOfMass ≈ referenceInertia.centerOfMass
            @test inertia.inertiaMatrix ≈ referenceInertia.inertiaMatrix
            @test length(pazy.beamNodes) == ne + 1
            @test all(iszero,pazy.model.gravityVector)
            Λ = pazy.parameters.sweepAngle
            weights = Matrix{Float64}(I,ne+1,ne+1)
            grid,positions,offsets = wing_geometry(pazy.model,c,spar,4,weights; sweepAngle=Λ)
            @test size(grid) == (3,5,ne+1)
            @test grid[:,end,end] - grid[:,1,end] ≈ [c,0.0,0.0] atol=1e-14
            # AeroBeams obtains nodal positions through numerical integration.
            @test grid[:,1,end] - grid[:,1,1] ≈ [L*sin(Λ),L*cos(Λ),0.0] atol=1e-9
            @test maximum(abs,grid[2,:,1]) < 1e-14
            @test maximum(abs,grid[3,:,:]) < 1e-14
            @test wingtip_twist_degrees(pazy.model) ≈ 0.0 atol=1e-12

            # UVLM nodal forces and mapped wrenches must have identical resultants.
            _,_,panels = UVLM.grid_to_surface_panels(grid)
            points = UVLM.imperial_nodal_positions(panels)
            forces = [[sin(k+j),cos(k-j),1.0] for k in axes(points,1),j in axes(points,2)]
            loads = beam_loads((forces=[forces],positions=[points]),offsets,weights)
            @test size(loads) == (6,ne+1)
            @test vec(sum(loads[1:3,:]; dims=2)) ≈ sum(A_TO_UVLM' * f for f in forces)
            aeroMoment = sum(cross(A_TO_UVLM' * points[i],A_TO_UVLM' * forces[i]) for i in eachindex(points))
            beamMoment = sum(cross(positions[:,n],loads[1:3,n]) + loads[4:6,n] for n in axes(positions,2))
            @test aeroMoment ≈ beamMoment
            apply_pazy_nodal_loads!((model=pazy.model,),loads,0.0)
            @test all(pazy.model.BCs[n+1].state.currentValue ≈ loads[:,n] for n in axes(loads,2))
        end
    end
    @test_throws AssertionError create_pazy_uvlm_structure(sweepAngle=deg2rad(5))
    @test_throws AssertionError create_pazy_uvlm_structure(sweepAngle=deg2rad(10),useAePW3Matrices=true)
    @test_throws AssertionError create_pazy_uvlm_structure(mesh=:refined,useAePW3Matrices=true)
    @test_throws AssertionError create_pazy_uvlm_structure(sweepAngle=deg2rad(10),withWingtipRod=false)
    @test_throws AssertionError create_pazy_uvlm_structure(withWingtipRod=false,beamWingtipRod=true)
    @test_throws AssertionError create_pazy_uvlm_structure(beamWingtipRod=true,sweptTipMassInertiaCorrection=true)
    @test_throws AssertionError create_pazy_uvlm_structure(sweepAngle=deg2rad(20),beamWingtipRod=true)
end

if "--dynamic" in ARGS
    options = (; airspeed=15.0,density=1.225,angle_of_attack=deg2rad(1.0),sideslip=0.0,
        initial_airspeed_fraction=1.0,airspeed_ramp_duration=0.0,
        duration=0.005,settling_time=0.0,perturbation_amplitude=0.01,
        perturbation_duration=0.01,chordwise_panels=2,spanwise_panels=3,
        symmetric_wing=true,maximum_wake_rows=4,core_radius=1e-3,
        wake_shedding_fraction=0.1,aerodynamic_interaction=true,
        save_uvlm_history=true,uvlm_save_frequency=1,
        newton_maximum_iterations=50,newton_absolute_tolerance=1e-8,
        newton_relative_tolerance=1e-8,newton_display_iterations=false,
        newton_always_update_jacobian=true,animation_frames=2,progress_frequency=100,
        coupling_maximum_iterations=40,coupling_relaxation=0.5)
    @testset "Coupled Pazy configurations" begin
        lastResult = nothing
        for pazyOptions in ((;useAePW3Matrices=true),(;withWingtipRod=false,mesh=:refined),
            (;sweepAngle=deg2rad(20)),(;sweepAngle=deg2rad(10),beamWingtipRod=true,sweptTipMassConfig=:LE))
            result = run_pazy_wing_uvlm(; options...,pazyOptions...)
            @test result.aerodynamic_problem.results.completed
            @test all(isfinite,result.tip_out_of_plane)
            @test all(isfinite,result.tip_twist_degrees)
            @test maximum(abs,result.tip_out_of_plane) > 0
            wing = result.structural_problem.model.beams[1]
            @test last(result.tip_out_of_plane) == -wing.elements[end].state.nodal.u2[1]
            tipRotation = AeroBeams.rotation_tensor_WM(wing.elements[end].state.nodal.p2_b)
            @test last(result.tip_twist_degrees) ≈ asind((tipRotation * [0.0,1.0,0.0])[3])
            @test size(result.aerodynamic_nodal_load_history,2) == wing.discretization.nElements+1
            @test maximum(result.coupling_geometry_residual) <= 1e-5
            @test maximum(result.coupling_load_residual) <= 1e-3
            @test result.aerodynamic_problem.model.reference.b ≈
                  wing.geometry.ℓ*cos(result.pazy_parameters.sweepAngle)
            lastResult = result
        end
        staticOptions = (; airspeed=options.airspeed,density=options.density,angle_of_attack=0.0,
            chordwise_panels=2,spanwise_panels=3,symmetric_wing=true,maximum_wake_rows=4,
            core_radius=1e-3,wake_shedding_fraction=0.1,aerodynamic_interaction=true,
            maximum_aerodynamic_time=0.1,minimum_aerodynamic_time=0.0,
            aerodynamic_window_time=0.01,aerodynamic_force_tolerance=1e-4,
            aerodynamic_circulation_tolerance=1e-4,maximum_static_iterations=5,
            static_relaxation=0.5,geometry_tolerance=1e-5,load_tolerance=1e-4,
            consecutive_equilibrium_iterations=2,newton_maximum_iterations=50,
            newton_absolute_tolerance=1e-8,newton_relative_tolerance=1e-8,
            newton_display_iterations=false,progress_frequency=1000)
        for pazyOptions in ((;withWingtipRod=false,mesh=:refined),
            (;sweepAngle=deg2rad(20),sweptTipMassConfig=:TE),
            (;sweepAngle=deg2rad(10),beamWingtipRod=true,sweptTipMassConfig=:LE))
            equilibrium = run_pazy_wing_uvlm_static(; staticOptions...,pazyOptions...)
            @test equilibrium.converged
            @test norm(equilibrium.tip_displacement_A) < 1e-10
            originalState = copy(equilibrium.structural_problem.systemArrays.x)
            seeded = run_pazy_wing_uvlm(; options...,pazyOptions...,angle_of_attack=0.0,
                perturbation_amplitude=0.0,initial_equilibrium=equilibrium)
            @test seeded.initial_equilibrium_state.structural_state == originalState
            @test equilibrium.structural_problem.systemArrays.x == originalState
            @test maximum(abs,seeded.tip_out_of_plane) < 1e-10
            @test seeded.pazy_parameters == equilibrium.pazy_parameters
            @test_throws AssertionError run_pazy_wing_uvlm(; options...,pazyOptions...,
                initial_equilibrium=equilibrium,rootTorsionCorr=false)
        end
        if "--animations" in ARGS
            include(joinpath(@__DIR__,"..","PazyWingUVLMVisualization.jl"))
            directory = mktempdir(joinpath(REPOSITORY_ROOT,"tmp"); cleanup=false)
            serialize(joinpath(directory,"swept_rod_result.jls"),lastResult)
            println("Swept-rod animation previews: $directory")
            @test isfile(save_structural_animation(lastResult;
                output_path=joinpath(directory,"swept_rod_structure.gif"),fps=2))
            @test isfile(save_wake_animation(lastResult;
                output_path=joinpath(directory,"swept_rod_wake.gif"),fps=2))
        end
    end
end
