using DelimitedFiles
using WingPropellerUVLM

"""
Configuration for the rigid rectangular-wing free-wake demonstration.

The default temporal resolution convects the wake approximately one quarter
chord per step. Environment variables documented in `examples/README.md` may be
used for convergence studies without editing this file.
"""
Base.@kwdef struct RectangularWingCase
    span::Float64 = 7.5
    chord::Float64 = 1.0
    angleOfAttackDegrees::Float64 = 5.0
    freestreamSpeed::Float64 = 100.0
    airDensity::Float64 = 1.225
    spanwisePanels::Int = 30
    chordwisePanels::Int = 4
    stepsPerConvectiveTime::Int = 4
    convectiveTimes::Float64 = 15.0
    averagingConvectiveTimes::Float64 = 5.0
    vtkStride::Int = 5
    outputDirectory::String = joinpath(@__DIR__, "output", "rectangular_wing_free_wake")
end

readEnv(::Type{T},name,default) where {T} = parse(T,get(ENV,name,string(default)))

function caseFromEnvironment()
    defaults = RectangularWingCase()
    return RectangularWingCase(
        span=readEnv(Float64,"UVLM_RECT_SPAN",defaults.span),
        chord=readEnv(Float64,"UVLM_RECT_CHORD",defaults.chord),
        angleOfAttackDegrees=readEnv(Float64,"UVLM_RECT_ALPHA_DEG",defaults.angleOfAttackDegrees),
        freestreamSpeed=readEnv(Float64,"UVLM_RECT_VINF",defaults.freestreamSpeed),
        airDensity=readEnv(Float64,"UVLM_RECT_RHO",defaults.airDensity),
        spanwisePanels=readEnv(Int,"UVLM_RECT_NS",defaults.spanwisePanels),
        chordwisePanels=readEnv(Int,"UVLM_RECT_NC",defaults.chordwisePanels),
        stepsPerConvectiveTime=readEnv(Int,"UVLM_RECT_STEPS_PER_CHORD",defaults.stepsPerConvectiveTime),
        convectiveTimes=readEnv(Float64,"UVLM_RECT_CONVECTIVE_TIMES",defaults.convectiveTimes),
        averagingConvectiveTimes=readEnv(Float64,"UVLM_RECT_AVERAGING_TIMES",defaults.averagingConvectiveTimes),
        vtkStride=readEnv(Int,"UVLM_RECT_VTK_STRIDE",defaults.vtkStride),
        outputDirectory=get(ENV,"UVLM_RECT_OUTPUT_DIR",defaults.outputDirectory),
    )
end

function validate(case::RectangularWingCase)
    case.span > 0 || error("span must be positive")
    case.chord > 0 || error("chord must be positive")
    case.freestreamSpeed > 0 || error("freestream speed must be positive")
    case.airDensity > 0 || error("air density must be positive")
    case.spanwisePanels >= 2 || error("at least two spanwise panels are required")
    case.chordwisePanels >= 1 || error("at least one chordwise panel is required")
    case.stepsPerConvectiveTime >= 1 || error("steps per convective time must be positive")
    case.convectiveTimes > 0 || error("simulation duration must be positive")
    0 < case.averagingConvectiveTimes <= case.convectiveTimes ||
        error("averaging duration must be within the simulation duration")
    case.vtkStride >= 1 || error("VTK stride must be positive")
    return case
end

function initializeSystem(case::RectangularWingCase,maximumWakeRows::Int)
    halfSpan = case.span/2
    xle = [0.0,0.0]
    yle = [-halfSpan,halfSpan]
    zle = [0.0,0.0]
    chord = fill(case.chord,2)
    twist = zeros(2)
    dihedral = zeros(2)

    grid,ratios = wing_to_grid(
        xle,yle,zle,chord,twist,dihedral,
        case.spanwisePanels,case.chordwisePanels;
        spacing_s=Cosine(),
        spacing_c=Uniform(),
    )
    _,returnedRatios,surface = grid_to_surface_panels(grid;ratios)

    system = System([grid];nw=[maximumWakeRows])
    system.ratios[1] .= returnedRatios
    system.surfaces[1] .= surface
    system.previous_surfaces[1] .= surface
    system.symmetric .= false
    system.surface_id .= 1
    system.wake_finite_core .= true
    system.trailing_vortices .= false

    reference = Reference(
        case.span*case.chord,
        case.chord,
        case.span,
        [case.chord/4,0.0,0.0],
        case.freestreamSpeed,
        case.airDensity,
    )
    freestream = Freestream(
        case.freestreamSpeed,
        deg2rad(case.angleOfAttackDegrees),
        0.0,
        zeros(3),
    )
    system.reference[] = reference
    system.freestream[] = freestream
    system.Γ .= 0

    return system,grid,freestream
end

function meanRows(values,rows)
    return vec(sum(view(values,rows,:);dims=1))/length(rows)
end

function writeCsv(path,header,columns...)
    count = length(columns[1])
    all(length(column) == count for column in columns) ||
        error("CSV columns must have equal length")
    open(path,"w") do io
        println(io,join(header,','))
        for i in 1:count
            println(io,join((column[i] for column in columns),','))
        end
    end
    return path
end

function svgPolyline(x,y,xMap,yMap)
    return join(("$(round(xMap(xi);digits=2)),$(round(yMap(yi);digits=2))" for (xi,yi) in zip(x,y)),' ')
end

"""Write a dependency-free SVG line plot of mean lift per unit span."""
function writeLiftPlot(path,yCoordinates,meanLiftPerSpan,case)
    width,height = 960,600
    left,right,top,bottom = 105,35,70,90
    plotWidth = width-left-right
    plotHeight = height-top-bottom
    xMinimum,xMaximum = extrema(yCoordinates)
    rawMinimum,rawMaximum = extrema(meanLiftPerSpan)
    yMinimum = min(0.0,1.05*rawMinimum)
    yMaximum = max(0.0,1.10*rawMaximum)
    yMaximum > yMinimum || (yMaximum = yMinimum+1)
    xMap(x) = left+(x-xMinimum)/(xMaximum-xMinimum)*plotWidth
    yMap(y) = top+(yMaximum-y)/(yMaximum-yMinimum)*plotHeight
    points = svgPolyline(yCoordinates,meanLiftPerSpan,xMap,yMap)

    xTicks = range(xMinimum,xMaximum;length=7)
    yTicks = range(yMinimum,yMaximum;length=6)
    open(path,"w") do io
        println(io,"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"$width\" height=\"$height\" viewBox=\"0 0 $width $height\">")
        println(io,"<rect width=\"100%\" height=\"100%\" fill=\"white\"/>")
        println(io,"<text x=\"$(width/2)\" y=\"32\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"22\">Mean spanwise lift distribution</text>")
        println(io,"<text x=\"$(width/2)\" y=\"55\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"14\">Rectangular wing: b=$(case.span) m, c=$(case.chord) m, alpha=$(case.angleOfAttackDegrees) deg, V=$(case.freestreamSpeed) m/s</text>")
        for value in yTicks
            pixel = yMap(value)
            println(io,"<line x1=\"$left\" y1=\"$pixel\" x2=\"$(left+plotWidth)\" y2=\"$pixel\" stroke=\"#dddddd\" stroke-width=\"1\"/>")
            println(io,"<text x=\"$(left-12)\" y=\"$(pixel+5)\" text-anchor=\"end\" font-family=\"sans-serif\" font-size=\"13\">$(round(value;digits=1))</text>")
        end
        for value in xTicks
            pixel = xMap(value)
            println(io,"<line x1=\"$pixel\" y1=\"$top\" x2=\"$pixel\" y2=\"$(top+plotHeight)\" stroke=\"#eeeeee\" stroke-width=\"1\"/>")
            println(io,"<text x=\"$pixel\" y=\"$(top+plotHeight+25)\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"13\">$(round(value;digits=2))</text>")
        end
        println(io,"<rect x=\"$left\" y=\"$top\" width=\"$plotWidth\" height=\"$plotHeight\" fill=\"none\" stroke=\"black\" stroke-width=\"1.5\"/>")
        println(io,"<polyline points=\"$points\" fill=\"none\" stroke=\"#005f9e\" stroke-width=\"3\" stroke-linejoin=\"round\"/>")
        for (x,y) in zip(yCoordinates,meanLiftPerSpan)
            println(io,"<circle cx=\"$(xMap(x))\" cy=\"$(yMap(y))\" r=\"3\" fill=\"#005f9e\"/>")
        end
        println(io,"<text x=\"$(left+plotWidth/2)\" y=\"$(height-28)\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"16\">Spanwise coordinate, y (m)</text>")
        println(io,"<text x=\"27\" y=\"$(top+plotHeight/2)\" text-anchor=\"middle\" transform=\"rotate(-90 27 $(top+plotHeight/2))\" font-family=\"sans-serif\" font-size=\"16\">Mean lift per unit span (N/m)</text>")
        println(io,"</svg>")
    end
    return path
end

function saveResults(
    case,system,grid,time,clHistory,sectionClHistory,surfaceHistory,
    propertyHistory,wakeHistory,vtkTimeSteps,Δt,averageRows,
)
    mkpath(case.outputDirectory)
    liftingLineCoordinates,liftingLineChords = lifting_line_geometry([grid])
    yEdges = vec(liftingLineCoordinates[1][2,:])
    yCenters = (yEdges[1:end-1]+yEdges[2:end])/2
    panelWidths = diff(yEdges)
    localChord = (liftingLineChords[1][1:end-1]+liftingLineChords[1][2:end])/2
    meanSectionCl = meanRows(sectionClHistory,averageRows)

    dynamicPressure = 0.5*case.airDensity*case.freestreamSpeed^2
    meanLiftPerSpan = dynamicPressure.*localChord.*meanSectionCl
    liftHistory = dynamicPressure*(case.span*case.chord).*clHistory
    meanCl = sum(@view(clHistory[averageRows]))/length(averageRows)
    integratedMeanLift = sum(meanLiftPerSpan.*panelWidths)
    integratedMeanCl = integratedMeanLift/(dynamicPressure*case.span*case.chord)

    writeCsv(
        joinpath(case.outputDirectory,"lift_history.csv"),
        ("time_s","CL","lift_N"),
        time,clHistory,liftHistory,
    )
    writeCsv(
        joinpath(case.outputDirectory,"mean_lift_distribution.csv"),
        ("y_m","two_y_over_span","mean_section_cl","mean_lift_per_span_N_per_m"),
        yCenters,2 .* yCenters ./ case.span,meanSectionCl,meanLiftPerSpan,
    )
    writeLiftPlot(
        joinpath(case.outputDirectory,"mean_lift_distribution.svg"),
        yCenters,meanLiftPerSpan,case,
    )

    savedSteps = findall(vtkTimeSteps)
    vtkIntervals = diff(vcat(0,savedSteps)).*Δt
    vtkBase = joinpath(case.outputDirectory,"rectangular_wing_free_wake")
    metadata = Dict(
        "span_m"=>case.span,
        "chord_m"=>case.chord,
        "alpha_deg"=>case.angleOfAttackDegrees,
        "freestream_mps"=>case.freestreamSpeed,
        "air_density_kgpm3"=>case.airDensity,
    )
    write_vtk(
        vtkBase,surfaceHistory,propertyHistory,wakeHistory,vtkIntervals;
        symmetric=[false],metadata,
    )

    summaryPath = joinpath(case.outputDirectory,"summary.txt")
    open(summaryPath,"w") do io
        println(io,"Rigid rectangular-wing UVLM free-wake case")
        println(io,"span_m = $(case.span)")
        println(io,"chord_m = $(case.chord)")
        println(io,"angle_of_attack_deg = $(case.angleOfAttackDegrees)")
        println(io,"freestream_speed_mps = $(case.freestreamSpeed)")
        println(io,"air_density_kgpm3 = $(case.airDensity)")
        println(io,"spanwise_panels = $(case.spanwisePanels)")
        println(io,"chordwise_panels = $(case.chordwisePanels)")
        println(io,"time_step_s = $Δt")
        println(io,"time_steps = $(length(time))")
        println(io,"simulated_time_s = $(time[end])")
        println(io,"averaging_start_s = $(time[first(averageRows)])")
        println(io,"averaging_end_s = $(time[last(averageRows)])")
        println(io,"mean_CL_from_body_forces = $meanCl")
        println(io,"mean_CL_from_sectional_integration = $integratedMeanCl")
        println(io,"mean_lift_N = $integratedMeanLift")
        println(io,"maximum_active_wake_rows = $(maximum(system.nwake))")
        println(io,"vtk_frames = $(length(surfaceHistory))")
    end

    return (;
        meanCl,
        integratedMeanCl,
        integratedMeanLift,
        summaryPath,
        vtkBase,
    )
end

function run(case::RectangularWingCase=caseFromEnvironment())
    validate(case)
    convectiveTime = case.chord/case.freestreamSpeed
    Δt = convectiveTime/case.stepsPerConvectiveTime
    numberOfSteps = round(Int,case.convectiveTimes*case.stepsPerConvectiveTime)
    averageStepCount = round(Int,case.averagingConvectiveTimes*case.stepsPerConvectiveTime)
    averageRows = (numberOfSteps-averageStepCount+1):numberOfSteps
    maximumWakeRows = numberOfSteps

    system,grid,freestream = initializeSystem(case,maximumWakeRows)
    repeatedPoints = repeated_trailing_edge_points(system.surfaces)
    activeWakeRows = [0]
    maximumWakeRowsBySurface = [maximumWakeRows]
    time = collect(1:numberOfSteps).*Δt
    clHistory = zeros(numberOfSteps)
    sectionClHistory = zeros(numberOfSteps,case.spanwisePanels)

    surfaceHistory = Vector{Vector{Matrix{SurfacePanel{Float64}}}}()
    propertyHistory = Vector{Vector{Matrix{PanelProperties{Float64}}}}()
    wakeHistory = Vector{Vector{Matrix{WakePanel{Float64}}}}()
    vtkTimeSteps = falses(numberOfSteps)

    println("Running rigid rectangular-wing free-wake UVLM case")
    println("  lattice: $(case.chordwisePanels) chordwise x $(case.spanwisePanels) spanwise panels")
    println("  dt = $Δt s, steps = $numberOfSteps, wake rows = $maximumWakeRows")

    for step in 1:numberOfSteps
        # This is the accepted-state pattern used by the run_chang partitioned
        # drivers. A coupled structural iteration may call advance_uvlm_trial!
        # repeatedly from this snapshot; only the accepted call is retained.
        copy_surfaces_to_previous!(system,1)
        snapshot = snapshot_uvlm(system)
        advance_uvlm_trial!(
            system,snapshot,freestream,Δt;
            repeatedPoints,
            activeWakeRows,
            η=0.1,
            calculateInfluenceMatrix=true,
            nearFieldAnalysis=true,
            derivatives=false,
        )
        commit_wake_rows!(activeWakeRows,maximumWakeRowsBySurface)
        system.nwake .= activeWakeRows

        totalForceCoefficient,_ = body_forces(system;frame=Wind())
        sectionalForceCoefficient,_ = lifting_line_coefficients(system;frame=Wind())
        clHistory[step] = totalForceCoefficient[3]
        sectionClHistory[step,:] .= sectionalForceCoefficient[1][3,:]

        saveVtk = step == 1 || step % case.vtkStride == 0 || step == numberOfSteps
        if saveVtk
            vtkTimeSteps[step] = true
            push!(surfaceHistory,[copy(surface) for surface in system.surfaces])
            push!(propertyHistory,[copy(properties) for properties in system.properties])
            push!(wakeHistory,[copy(system.wakes[i][1:activeWakeRows[i],:]) for i in eachindex(system.wakes)])
        end

        if step == 1 || step % max(1,div(numberOfSteps,10)) == 0 || step == numberOfSteps
            println("  step $step/$numberOfSteps: t=$(round(time[step];digits=5)) s, CL=$(round(clHistory[step];digits=6)), wake rows=$(activeWakeRows[1])")
        end
    end

    results = saveResults(
        case,system,grid,time,clHistory,sectionClHistory,surfaceHistory,
        propertyHistory,wakeHistory,vtkTimeSteps,Δt,averageRows,
    )
    println("Completed: mean CL = $(results.meanCl)")
    println("Integrated sectional mean CL = $(results.integratedMeanCl)")
    println("Mean lift = $(results.integratedMeanLift) N")
    println("Outputs: $(abspath(case.outputDirectory))")
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end
