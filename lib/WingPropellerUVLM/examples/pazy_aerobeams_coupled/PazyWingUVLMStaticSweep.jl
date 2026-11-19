# Spanwise normal forces and pitching moments about nonlinear static Pazy equilibria.
if !isdefined(@__MODULE__,:run_pazy_wing_uvlm_static)
    include(joinpath(@__DIR__,"PazyWingUVLMStatic.jl"))
end
using DelimitedFiles

"""Sectional cn and cm at spanwise nodes of the explicitly modelled semispan.

cn = dot(F_node,n_node)/(0.5*rho*U^2*c*dy0), where dy0 is the
UNDEFORMED projected tributary width: one panel width at interior nodes,
half at the root and tip. F_node sums chordwise vortex-vertex forces after
the same conservative redistribution used for structural load transfer.
n_node follows the physical chord and the centered span tangent (one-sided
at boundaries), pointing upward on the undeformed wing. Endpoint cn values
are half-strip averages, not pointwise tip/root limits. This is a finite-wing
coefficient referenced to freestream dynamic pressure, including induced-flow effects.

cm is the nodal wrench about the deformed physical quarter chord or spar
divided by q_inf * chord^2 * dy0. Moments use the UVLM vortex vertices.
The pitch axis is n_node × t_chord (nose-up positive). Supply spar_fraction
for the native spar reference, measured from the physical leading edge.
alpha_local is the signed freestream incidence in the deformed chord/normal
plane. It includes elastic deformation, but excludes induced downwash from
the angle definition. cn_alpha_local = cn/alpha_local is a secant ratio, not
a tangent derivative. Angles too close to zero give NaN ratios.
"""
function pazy_spanwise_cn(problem,span,sweepAngle; local_alpha_tolerance=1e-8,spar_fraction=nothing)
    system = problem.state.system
    grid = system.grids[1]
    nc,ns = size(system.surfaces[1])
    reference = problem.model.reference
    eta = collect(0:ns) ./ ns
    tributaryArea = fill(reference.c * span * cos(sweepAngle) / ns,ns + 1)
    tributaryArea[1] /= 2
    tributaryArea[end] /= 2
    sectionForces = zeros(3,ns + 1)
    sectionNormals = zeros(3,ns + 1)
    sectionPositions = zeros(3,ns + 1)
    cn = zeros(ns + 1)
    cmQuarterChord = zeros(ns + 1)
    cmSpar = fill(NaN,ns + 1)
    alphaLocal = fill(NaN,ns + 1)
    localRatio = fill(NaN,ns + 1)
    localAlphaValid = falses(ns + 1)
    velocity = UVLM.freestream_velocity(system.freestream[])
    q = reference.rho * reference.V^2 / 2
    # Segment forces split between their endpoints; panel unsteady forces
    # split between four corners. Sum those dimensional nodal forces along
    # each chord before projecting and dividing by its tributary area.
    nodalForces = UVLM.imperial_nodal_forces(system.span_seg_forces[1],
        system.chord_seg_forces[1],system.unsteady_forces[1])
    vortexPositions = UVLM.imperial_nodal_positions(system.surfaces[1])
    for j = 1:ns+1
        leadingEdge = UVLM.SVector{3}(grid[:,1,j])
        chordDirection = UVLM.SVector{3}(grid[:,end,j]) - leadingEdge
        quarterChord = leadingEdge + chordDirection / 4
        moment = zero(quarterChord)
        for i = 1:nc+1
            force = nodalForces[i,j]
            sectionForces[:,j] .+= force
            moment += cross(vortexPositions[i,j] - quarterChord,force)
        end
        left,right = max(1,j - 1),min(ns + 1,j + 1)
        spanDirection = (grid[:,1,right] + grid[:,end,right] - grid[:,1,left] - grid[:,end,left]) / 2
        normalDirection = cross(chordDirection,spanDirection)
        norm(normalDirection) > eps(Float64) || error("Degenerate aerodynamic section $j")
        normal = normalize(normalDirection)
        chordUnit = normalize(chordDirection)
        pitchAxis = cross(normal,chordUnit)
        sectionNormals[:,j] = normal
        sectionPositions[:,j] = leadingEdge + chordDirection / 2
        cn[j] = dot(sectionForces[:,j],normal) / (q * tributaryArea[j])
        cmQuarterChord[j] = dot(moment,pitchAxis) / (q * tributaryArea[j] * reference.c)
        if !isnothing(spar_fraction)
            # Shift the moment along the deformed physical chord to the spar.
            cmSpar[j] = cmQuarterChord[j] +
                (spar_fraction - 0.25) * (norm(chordDirection) / reference.c) * cn[j]
        end
        chordVelocity = dot(velocity,chordUnit)
        normalVelocity = dot(velocity,normal)
        if hypot(chordVelocity,normalVelocity) > eps(Float64) * norm(velocity)
            alphaLocal[j] = atan(normalVelocity,chordVelocity)
            localAlphaValid[j] = abs(alphaLocal[j]) > local_alpha_tolerance
            localAlphaValid[j] && (localRatio[j] = cn[j] / alphaLocal[j])
        end
    end
    return (; eta,s_m=eta .* span,y_reference_m=eta .* (span * cos(sweepAngle)),
        cn,cm_quarter_chord=cmQuarterChord,cm_spar=cmSpar,
        strip_area_m2=tributaryArea,forces_UVLM=sectionForces,
        normals_UVLM=sectionNormals,positions_UVLM=sectionPositions,
        alpha_local=alphaLocal,cn_alpha_local=localRatio,local_alpha_valid=localAlphaValid)
end

function pazy_static_derivatives(sections,alphaStep; local_alpha_tolerance=1e-8)
    minus,plus = sections
    deltaAlpha = 2 * alphaStep
    deltaCn = plus.cn - minus.cn
    deltaAlphaLocal = rem2pi.(plus.alpha_local - minus.alpha_local,RoundNearest)
    cnAlphaLocal = fill(NaN,length(deltaCn))
    valid = falses(length(deltaCn))
    for j in eachindex(deltaCn)
        if isfinite(deltaAlphaLocal[j]) && abs(deltaAlphaLocal[j]) > local_alpha_tolerance
            slope = deltaCn[j] / deltaAlphaLocal[j]
            valid[j] = isfinite(slope)
            valid[j] && (cnAlphaLocal[j] = slope)
        end
    end
    localIncidence = (; delta_cn=deltaCn,delta_alpha_rad=deltaAlphaLocal,cn_alpha=cnAlphaLocal,valid)
    return (; cn=deltaCn / deltaAlpha,
        cm_quarter_chord=(plus.cm_quarter_chord - minus.cm_quarter_chord) / deltaAlpha,
        cm_spar=(plus.cm_spar - minus.cm_spar) / deltaAlpha,
        local_incidence=localIncidence)
end

# Reconstruct the existing static solver's material interpolation, retaining
# the native spar reference and streamwise chord orientation for swept wings.
function pazy_static_interface(equilibrium)
    parameters = equilibrium.pazy_parameters
    ne,span,chord,spar,beamNodes = AeroBeams.geometrical_properties_Pazy(;
        mesh=parameters.mesh,Λ=parameters.sweepAngle,useAePW3Matrices=parameters.useAePW3Matrices)
    problem = equilibrium.aerodynamic_problem
    nc,ns = size(problem.state.system.surfaces[1])
    weights = zeros(ns + 1,ne + 1)
    for j = 1:ns+1
        eta = (j - 1) / ns
        left = min(searchsortedlast(beamNodes,eta),ne)
        fraction = (eta - beamNodes[left]) / (beamNodes[left+1] - beamNodes[left])
        weights[j,left],weights[j,left+1] = 1 - fraction,fraction
    end
    grid,_,offsets = wing_geometry(equilibrium.structural_problem.model,chord,spar,nc,weights;
        sweepAngle=parameters.sweepAngle)
    return (; grid,offsets,weights,span,chord)
end

"""Sweep static equilibria and obtain spanwise cn, cn/alpha_local and cm.

The default :local_ratio uses one static equilibrium per condition, with
prescribed incidence held constant. alpha_local follows the deformed physical
chord and surface normal; it excludes induced downwash. cn/alpha_local is a
secant approximation to a lift slope and is NaN near zero local incidence.
Optional derivative modes also give both moment derivatives. Angles and alpha_step are radians;
these derivatives are per radian of imposed freestream
incidence. :equilibrium re-solves the flexible structure at alpha +/- step;
:fixed_shape holds the central equilibrium geometry fixed at both incidences.
Both modes develop fresh free wakes. Neither is a 2D airfoil slope with
respect to local effective incidence. The same perturbations also give
delta cn / delta alpha_local for each station, measured along the imposed
global incidence change. static_options supplies the existing
run_pazy_wing_uvlm_static controls, except speed/incidence/sweep and wake times.
Wake times are specified here in convective units U*t/c, consistently at all
speeds. Unconverged cases retain diagnostics but export NaN slopes.

With warm_start=true each equilibrium starts from the previous one at the same sweep and
incidence angles, if the dynamic pressure ratio is within [1/5,5] (structure and applied loads
scaled by that ratio; otherwise the case starts cold),
and the perturbed solves of a case start from its central equilibrium (and reuse
its wake if static_options has wake_reuse=true). The converged equilibria satisfy
the same tolerances; only the path to them changes. See run_pazy_wing_uvlm_static
for relaxation_mode, geometry_cap_chords and wake_reuse.
"""
function run_pazy_wing_uvlm_static_sweep(; airspeeds,angles_of_attack,
    sweep_angles=[0.0],alpha_step=deg2rad(0.25),
    slope_modes=(:local_ratio,),static_options,local_alpha_tolerance=1e-8,
    maximum_aerodynamic_chords=120.0,minimum_aerodynamic_chords=25.0,
    aerodynamic_window_chords=4.0,output_directory=nothing,save_plots=true,warm_start=true)

    for (name,values) in (("airspeeds",airspeeds),("angles_of_attack",angles_of_attack),
        ("sweep_angles",sweep_angles))
        !isempty(values) && all(isfinite,values) || throw(ArgumentError("$name must be nonempty and finite"))
        length(unique(values)) == length(values) || throw(ArgumentError("$name must contain distinct values"))
    end
    all(>(0),airspeeds) || throw(ArgumentError("UVLM requires strictly positive airspeeds"))
    all(angle -> abs(angle) < π / 2,angles_of_attack) ||
        throw(ArgumentError("Prescribed incidence must be between -π/2 and π/2 radians for forward flow"))
    !isempty(slope_modes) && all(mode -> mode in (:local_ratio,:equilibrium,:fixed_shape),slope_modes) ||
        throw(ArgumentError("Choose :local_ratio, :equilibrium and/or :fixed_shape"))
    derivativeRequested = any(mode -> mode in (:equilibrium,:fixed_shape),slope_modes)
    if derivativeRequested
        isfinite(alpha_step) && alpha_step > 0 || throw(ArgumentError("alpha_step must be positive and finite for derivatives"))
        all(angle -> abs(angle) + alpha_step < π / 2,angles_of_attack) ||
            throw(ArgumentError("Perturbed incidences must also be between -π/2 and π/2 radians"))
    end
    isfinite(local_alpha_tolerance) && local_alpha_tolerance > 0 ||
        throw(ArgumentError("local_alpha_tolerance must be positive and finite"))
    all(isfinite,(maximum_aerodynamic_chords,minimum_aerodynamic_chords,aerodynamic_window_chords)) &&
        maximum_aerodynamic_chords > 0 && minimum_aerodynamic_chords >= 0 && aerodynamic_window_chords > 0 ||
        throw(ArgumentError("Invalid convective wake times"))
    reserved = (:airspeed,:angle_of_attack,:sweepAngle,:maximum_aerodynamic_time,
        :minimum_aerodynamic_time,:aerodynamic_window_time)
    any(key -> haskey(static_options,key),reserved) &&
        throw(ArgumentError("Specify speed, incidence, sweep and wake times through the sweep keywords"))

    cases = NamedTuple[]
    for sweepAngle in sweep_angles
        # Validate native data availability before starting this sweep.
        _,span,chord,sparFraction,_ = AeroBeams.geometrical_properties_Pazy(; Λ=sweepAngle,
            mesh=get(static_options,:mesh,:nominal),useAePW3Matrices=get(static_options,:useAePW3Matrices,false))
        previous = nothing
        for angle_of_attack in angles_of_attack,airspeed in airspeeds
            println("\nStatic spanwise slope: sweep=$(rad2deg(sweepAngle)) deg, " *
                "alpha=$(rad2deg(angle_of_attack)) deg, U=$airspeed m/s")
            options = (; static_options...,airspeed,angle_of_attack,sweepAngle,
                maximum_aerodynamic_time=maximum_aerodynamic_chords * chord / airspeed,
                minimum_aerodynamic_time=minimum_aerodynamic_chords * chord / airspeed,
                aerodynamic_window_time=aerodynamic_window_chords * chord / airspeed)
            startArguments = (;)
            # Only the previous equilibrium at the SAME prescribed incidence and a comparable dynamic pressure
            # (ratio within [1/5,5]) is a good starting point: otherwise (e.g. the first speed of a new incidence,
            # which would inherit the largest deformation of the previous incidence, or a jump from 0.01 m/s)
            # the case starts cold
            if warm_start && !isnothing(previous) && previous.angle_of_attack == angle_of_attack &&
                0.2 <= (airspeed / previous.airspeed)^2 <= 5
                # Dynamic-pressure scaling as predictor; the wake cannot be reused (time step changes with speed)
                startArguments = (; warm_start=(; previous.warm_start...,aero=nothing),
                    warm_start_target_loads=previous.warm_start.loads .* (airspeed / previous.airspeed)^2)
            end
            equilibrium = run_pazy_wing_uvlm_static(; options...,startArguments...)
            # Interface of the central shape, taken before any perturbed solve modifies the structure
            interface = :fixed_shape in slope_modes && equilibrium.converged ? pazy_static_interface(equilibrium) : nothing
            # Perturbed solves chain from the central state: same structure, loads and (if wake_reuse) wake
            perturbationStart = equilibrium.warm_start
            section = pazy_spanwise_cn(equilibrium.aerodynamic_problem,span,sweepAngle;
                local_alpha_tolerance,spar_fraction=sparFraction)
            localRatioValid = section.local_alpha_valid .& (equilibrium.converged && :local_ratio in slope_modes)
            localRatio = [localRatioValid[j] ? section.cn_alpha_local[j] : NaN for j in eachindex(section.cn)]
            noSlope = fill(NaN,length(section.cn))
            equilibriumSlope = (; cn=copy(noSlope),cm_quarter_chord=copy(noSlope),cm_spar=copy(noSlope),
                local_incidence=nothing)
            fixedShapeSlope = (; cn=copy(noSlope),cm_quarter_chord=copy(noSlope),cm_spar=copy(noSlope),
                local_incidence=nothing)
            diagnostics = NamedTuple[]
            push!(diagnostics,(; solve=:central,converged=equilibrium.converged,
                reason=equilibrium.termination_reason,history=equilibrium.history,
                history_columns=equilibrium.history_columns))
            equilibriumConverged = fixedShapeConverged = false
            if equilibrium.converged
                if :equilibrium in slope_modes
                    sections = typeof(section)[]
                    flags = Bool[]
                    for (label,sign) in ((:equilibrium_minus,-1),(:equilibrium_plus,1))
                        perturbed = run_pazy_wing_uvlm_static(; options...,
                            angle_of_attack=angle_of_attack + sign * alpha_step,
                            (warm_start ? (; warm_start=perturbationStart) : (;))...)
                        perturbationStart = perturbed.warm_start
                        push!(sections,pazy_spanwise_cn(perturbed.aerodynamic_problem,span,sweepAngle;
                            spar_fraction=sparFraction))
                        push!(flags,perturbed.converged)
                        push!(diagnostics,(; solve=label,converged=perturbed.converged,
                            reason=perturbed.termination_reason,history=perturbed.history,
                            history_columns=perturbed.history_columns))
                    end
                    equilibriumConverged = all(flags)
                    equilibriumConverged &&
                        (equilibriumSlope = pazy_static_derivatives(sections,alpha_step; local_alpha_tolerance))
                end
                if :fixed_shape in slope_modes
                    problem = equilibrium.aerodynamic_problem
                    # Fixed shape: the central wake is a good initial condition for the perturbed incidences
                    fixedShapeAero = warm_start && get(options,:wake_reuse,false) ? perturbationStart.aero : nothing
                    sections = typeof(section)[]
                    flags = Bool[]
                    for (label,sign) in ((:fixed_shape_minus,-1),(:fixed_shape_plus,1))
                        point = UVLM.create_OperatingPoint(; airspeed,density=options.density,
                            angleOfAttack=angle_of_attack + sign * alpha_step,sideslip=0.0)
                        wake = settle_static_uvlm(interface.grid,interface.offsets,interface.weights,
                            problem.model.reference,point,problem.aeroSolver;
                            symmetric_wing=options.symmetric_wing,maximum_wake_rows=options.maximum_wake_rows,
                            dt=equilibrium.dt,maximum_time=options.maximum_aerodynamic_time,
                            minimum_time=options.minimum_aerodynamic_time,window_time=options.aerodynamic_window_time,
                            force_tolerance=options.aerodynamic_force_tolerance,
                            circulation_tolerance=options.aerodynamic_circulation_tolerance,
                            progress_frequency=options.progress_frequency,aero=fixedShapeAero)
                        fixedShapeAero = warm_start && get(options,:wake_reuse,false) ? wake.aerodynamic_problem : nothing
                        push!(sections,pazy_spanwise_cn(wake.aerodynamic_problem,span,sweepAngle;
                            spar_fraction=sparFraction))
                        push!(flags,wake.converged)
                        push!(diagnostics,(; solve=label,converged=wake.converged,
                            reason=wake.converged ? "Fixed-shape wake settled" : "Fixed-shape wake did not settle",
                            wake_time_s=wake.time,force_residual=wake.force_residual,
                            circulation_residual=wake.circulation_residual))
                    end
                    fixedShapeConverged = all(flags)
                    fixedShapeConverged &&
                        (fixedShapeSlope = pazy_static_derivatives(sections,alpha_step; local_alpha_tolerance))
                end
            end
            # Keep compact sectional data and convergence histories, rather than
            # every large structural/wake problem across the parameter sweep.
            case = (; airspeed,angle_of_attack,sweepAngle,
                alpha_step=derivativeRequested ? alpha_step : NaN,section,
                cn_alpha_local=localRatio,cn_alpha_local_valid=localRatioValid,
                cn_alpha_equilibrium=equilibriumSlope.cn,cn_alpha_fixed_shape=fixedShapeSlope.cn,
                cn_alpha_local_perturbations=(; equilibrium=equilibriumSlope.local_incidence,
                    fixed_shape=fixedShapeSlope.local_incidence),
                cm_quarter_chord_alpha_equilibrium=equilibriumSlope.cm_quarter_chord,
                cm_spar_alpha_equilibrium=equilibriumSlope.cm_spar,
                cm_quarter_chord_alpha_fixed_shape=fixedShapeSlope.cm_quarter_chord,
                cm_spar_alpha_fixed_shape=fixedShapeSlope.cm_spar,
                central_converged=equilibrium.converged,
                equilibrium_slope_converged=equilibriumConverged,
                fixed_shape_slope_converged=fixedShapeConverged,
                tip_out_of_plane=equilibrium.tip_out_of_plane,tip_twist_degrees=equilibrium.tip_twist_degrees,
                pazy_parameters=equilibrium.pazy_parameters,diagnostics)
            push!(cases,case)
            previous = equilibrium.converged ?
                (; airspeed,angle_of_attack,warm_start=(; equilibrium.warm_start...,aero=nothing)) : nothing
            if !all(d.converged for d in diagnostics)
                @warn "Unconverged static/wake solves: inspect status.csv; corresponding slopes are NaN" airspeed sweepAngle angle_of_attack
            end
            if !isnothing(output_directory)
                save_pazy_static_sweep((; cases,slope_modes,static_options,local_alpha_tolerance,
                    maximum_aerodynamic_chords,minimum_aerodynamic_chords,aerodynamic_window_chords);
                    output_directory,save_plots=false)
            end
        end
    end
    result = (; cases,slope_modes,static_options,local_alpha_tolerance,maximum_aerodynamic_chords,
        minimum_aerodynamic_chords,aerodynamic_window_chords)
    isnothing(output_directory) || save_pazy_static_sweep(result; output_directory,save_plots)
    return result
end

function pazy_static_angle_folder(angle)
    degrees = rad2deg(angle)
    roundedDegrees = round(degrees; digits=10)
    # Remove degree/radian roundoff for ordinary input angles, while retaining
    # precision for angles that were not specified by that rounded value.
    displayDegrees = deg2rad(roundedDegrees) == angle ? roundedDegrees : degrees
    number = replace(string(displayDegrees),r"\.0$" => "")
    return "aoa$(replace(number,"." => "p","-" => "m"))deg"
end

function save_pazy_static_sweep(result; output_directory,save_plots=true)
    for angle in unique([case.angle_of_attack for case in result.cases])
        selected = filter(case -> case.angle_of_attack == angle,result.cases)
        angleResult = merge(result,(; cases=selected))
        save_pazy_static_sweep_cases(angleResult;
            output_directory=joinpath(output_directory,pazy_static_angle_folder(angle)),save_plots)
    end
    return output_directory
end

# Each angle directory contains its own tables, case numbering and figures.
function save_pazy_static_sweep_cases(result; output_directory,save_plots=true)
    mkpath(output_directory)
    rows = Vector{Any}[]
    momentRows = Vector{Any}[]
    perturbationRows = Vector{Any}[]
    statusRows = Vector{Any}[]
    derivativeModes = filter(mode -> mode in result.slope_modes,(:equilibrium,:fixed_shape))
    for (index,case) in enumerate(result.cases)
        section = case.section
        perturbations = get(case,:cn_alpha_local_perturbations,(;))
        for j in eachindex(section.cn)
            push!(rows,Any[index,case.airspeed,rad2deg(case.angle_of_attack),rad2deg(case.sweepAngle),
                rad2deg(case.alpha_step),section.eta[j],section.s_m[j],section.y_reference_m[j],
                section.positions_UVLM[:,j]...,case.central_converged ? section.cn[j] : NaN,
                case.cn_alpha_equilibrium[j],case.cn_alpha_fixed_shape[j],case.central_converged,
                case.equilibrium_slope_converged,case.fixed_shape_slope_converged,
                case.tip_out_of_plane,case.tip_twist_degrees,rad2deg(section.alpha_local[j]),
                section.alpha_local[j],case.cn_alpha_local[j],case.cn_alpha_local_valid[j]])
            push!(momentRows,Any[index,case.airspeed,rad2deg(case.angle_of_attack),rad2deg(case.sweepAngle),
                rad2deg(case.alpha_step),section.eta[j],section.s_m[j],section.y_reference_m[j],
                case.central_converged ? section.cm_quarter_chord[j] : NaN,
                case.central_converged ? section.cm_spar[j] : NaN,
                case.cm_quarter_chord_alpha_equilibrium[j],case.cm_spar_alpha_equilibrium[j],
                case.cm_quarter_chord_alpha_fixed_shape[j],case.cm_spar_alpha_fixed_shape[j],
                case.central_converged,case.equilibrium_slope_converged,case.fixed_shape_slope_converged])
            for mode in derivativeModes
                perturbation = get(perturbations,mode,nothing)
                values = isnothing(perturbation) ? (NaN,NaN,NaN,false) :
                    (perturbation.delta_cn[j],perturbation.delta_alpha_rad[j],
                     perturbation.cn_alpha[j],perturbation.valid[j])
                push!(perturbationRows,Any[index,case.airspeed,rad2deg(case.angle_of_attack),rad2deg(case.sweepAngle),
                    rad2deg(case.alpha_step),section.eta[j],section.s_m[j],section.y_reference_m[j],
                    mode,values...])
            end
        end
        for diagnostic in case.diagnostics
            push!(statusRows,Any[index,case.airspeed,rad2deg(case.angle_of_attack),rad2deg(case.sweepAngle),
                diagnostic.solve,diagnostic.converged,diagnostic.reason,
                get(diagnostic,:wake_time_s,NaN),get(diagnostic,:force_residual,NaN),
                get(diagnostic,:circulation_residual,NaN)])
            if haskey(diagnostic,:history)
                writedlm(joinpath(output_directory,"case_$(index)_$(diagnostic.solve)_convergence.csv"),
                    vcat(permutedims(diagnostic.history_columns),diagnostic.history),',')
            end
        end
    end
    columns = ["case","U_m_s","alpha_deg","sweep_deg","alpha_step_deg","eta","s_reference_m",
        "y_reference_m","x_deformed_UVLM_m","y_deformed_UVLM_m","z_deformed_UVLM_m","cn",
        "cn_alpha_equilibrium_per_rad","cn_alpha_fixed_shape_per_rad","central_converged",
        "equilibrium_slope_converged","fixed_shape_slope_converged","tip_OOP_m","tip_twist_deg",
        "alpha_local_deg","alpha_local_rad","cn_alpha_local_per_rad","cn_alpha_local_valid"]
    writedlm(joinpath(output_directory,"spanwise_cn_alpha.csv"),
        vcat(permutedims(columns),permutedims(hcat(rows...))),',')
    momentColumns = [columns[1:8]...,"cm_quarter_chord","cm_spar",
        "cm_quarter_chord_alpha_equilibrium_per_rad","cm_spar_alpha_equilibrium_per_rad",
        "cm_quarter_chord_alpha_fixed_shape_per_rad","cm_spar_alpha_fixed_shape_per_rad",
        "central_converged","equilibrium_slope_converged","fixed_shape_slope_converged"]
    writedlm(joinpath(output_directory,"spanwise_cm.csv"),
        vcat(permutedims(momentColumns),permutedims(hcat(momentRows...))),',')
    perturbationColumns = [columns[1:8]...,"mode","delta_cn","delta_alpha_local_rad",
        "cn_alpha_local_perturbation_per_rad","valid"]
    perturbationTable = permutedims(perturbationColumns)
    if !isempty(perturbationRows)
        perturbationTable = vcat(perturbationTable,permutedims(hcat(perturbationRows...)))
    end
    writedlm(joinpath(output_directory,"spanwise_cn_alpha_perturbations.csv"),perturbationTable,',')
    statusColumns = ["case","U_m_s","alpha_deg","sweep_deg","solve","converged","reason",
        "wake_time_s","wake_force_residual","wake_circulation_residual"]
    writedlm(joinpath(output_directory,"status.csv"),
        vcat(permutedims(statusColumns),permutedims(hcat(statusRows...))),',')
    open(joinpath(output_directory,"settings.txt"),"w") do io
        println(io,"Prescribed angle of attack [deg]: ",rad2deg(first(result.cases).angle_of_attack))
        println(io,"Default :local_ratio uses cn/alpha_local at one static equilibrium, without perturbations.")
        println(io,"alpha_local = atan(V_inf dot n_deformed, V_inf dot chord_downstream); radians.")
        println(io,"Local incidence includes deformation; induced downwash is excluded from its definition.")
        println(io,"Ratios are secant approximations, not tangent derivatives; units 1/rad.")
        println(io,"Local incidence cutoff [rad]: ",result.local_alpha_tolerance)
        println(io,"Optional derivative modes use central differences at alpha +/- alpha_step.")
        if !isempty(derivativeModes)
            println(io,"The perturbation CSV divides delta cn by the local incidence change along the same global perturbation.")
            println(io,"Local differences at or below the incidence cutoff have valid=false.")
        end
        println(io,"cn = F_node dot n_deformed / (q_inf * chord * undeformed projected tributary width).")
        println(io,"Section forces sum transferred vortex-vertex forces along each chord; total force is conserved.")
        println(io,"cm = moment about the deformed quarter chord or native spar / (q_inf * tributary area * reference chord).")
        println(io,"Moment arms start at the reference point and end at the force-bearing UVLM vortex vertices.")
        println(io,"The pitch axis is n_deformed cross chord_downstream. Positive cm is nose-up.")
        println(io,"Samples are at spanwise nodes, eta=0:1/spanwise_panels:1, including root and tip.")
        println(io,"Root/tip tributary widths are half a panel; endpoint cn values are half-strip averages.")
        println(io,"Local normals use the physical nodal chord and centered span tangent, one-sided at root/tip.")
        println(io,"eta follows undeformed spar arclength; positive sweep is the native aft sweep.")
        println(io,"UVLM axes: x downstream, y spanwise, z up. Leading edge upstream; wake downstream.")
        println(io,"symmetric_wing=true includes the image; forces describe only the explicit semispan.")
        println(io,"Requested slope modes: ",result.slope_modes)
        println(io,"Static options: ",result.static_options)
        println(io,"Wake U*t/c (maximum, minimum, window): ",
            (result.maximum_aerodynamic_chords,result.minimum_aerodynamic_chords,result.aerodynamic_window_chords))
    end
    if save_plots
        plotter = AeroBeams.Plots
        for sweepAngle in unique([case.sweepAngle for case in result.cases]),
            angle in unique([case.angle_of_attack for case in result.cases])
            selected = filter(case -> case.sweepAngle == sweepAngle && case.angle_of_attack == angle,result.cases)
            p1 = plotter.plot(xlabel="Normalized span eta = s/L (0: root, 1: tip)",
                ylabel="cn",xlims=(0.0,1.0),xticks=0.0:0.1:1.0,legend=:best)
            ratioOnly = all(==(:local_ratio),result.slope_modes)
            p2 = plotter.plot(xlabel="Normalized span eta = s/L (0: root, 1: tip)",
                ylabel=ratioOnly ? "cn / alpha_local [1/rad]" : "cn_alpha [1/rad]",
                xlims=(0.0,1.0),xticks=0.0:0.1:1.0,legend=:best)
            pQuarterChord = plotter.plot(xlabel="Normalized span eta = s/L (0: root, 1: tip)",
                ylabel="cm about quarter chord",xlims=(0.0,1.0),xticks=0.0:0.1:1.0,legend=:best)
            pSpar = plotter.plot(xlabel="Normalized span eta = s/L (0: root, 1: tip)",
                ylabel="cm about spar",xlims=(0.0,1.0),xticks=0.0:0.1:1.0,legend=:best)
            for (i,case) in enumerate(selected)
                if case.central_converged
                    label = "U=$(case.airspeed) m/s"
                    plotter.plot!(p1,case.section.eta,case.section.cn; label,color=i,lw=2)
                    plotter.plot!(pQuarterChord,case.section.eta,case.section.cm_quarter_chord; label,color=i,lw=2)
                    plotter.plot!(pSpar,case.section.eta,case.section.cm_spar; label,color=i,lw=2)
                end
                if :local_ratio in result.slope_modes && case.central_converged
                    plotter.plot!(p2,case.section.eta,case.cn_alpha_local;
                        label="U=$(case.airspeed), cn/alpha_local",color=i,lw=2,
                        linestyle=ratioOnly ? :solid : :dot)
                end
                case.equilibrium_slope_converged && plotter.plot!(p2,case.section.eta,case.cn_alpha_equilibrium;
                    label="U=$(case.airspeed), equilibrium",color=i,lw=2)
                case.fixed_shape_slope_converged && plotter.plot!(p2,case.section.eta,case.cn_alpha_fixed_shape;
                    label="U=$(case.airspeed), fixed shape",color=i,lw=2,linestyle=:dash)
            end
            title = "Pazy static UVLM, sweep=$(round(rad2deg(sweepAngle); digits=3)) deg, alpha=$(round(rad2deg(angle); sigdigits=6)) deg"
            figure = plotter.plot(p1,p2; layout=(2,1),size=(950,800),plot_title=title)
            momentFigure = plotter.plot(pQuarterChord,pSpar; layout=(2,1),size=(950,800),plot_title=title)
            # Full floating-point strings avoid collisions for nearby angles.
            suffix = replace("sweep$(rad2deg(sweepAngle))_alpha$(rad2deg(angle))","."=>"p","-"=>"m")
            plotter.savefig(figure,joinpath(output_directory,"spanwise_$(suffix).png"))
            plotter.savefig(momentFigure,joinpath(output_directory,"spanwise_cm_$(suffix).png"))
        end
    end
    return output_directory
end
