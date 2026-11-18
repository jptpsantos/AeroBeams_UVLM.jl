ENV["GKSwstype"] = "100"
include(joinpath(@__DIR__, "..", "lib", "WingPropellerUVLM", "examples",
    "pazy_aerobeams_coupled", "PazyWingUVLMCoupling.jl"))
include(joinpath(@__DIR__, "..", "lib", "WingPropellerUVLM", "examples",
    "pazy_aerobeams_coupled", "PazyWingUVLMVisualization.jl"))
using Printf, Serialization

function render_equal_scale_wake()
    output = joinpath(@__DIR__, "s20_wake_equal_scale")
    mkpath(output)
    U, alpha, sweep = 45.0, deg2rad(3.0), deg2rad(20.0)
    nc, ns, steps, visible_rows = 10, 20, 60, 10
    pazy = create_pazy_uvlm_structure(; sweepAngle=sweep)
    (; model,n_elements,span,chord,sparFraction,beamNodes) = pazy
    weights = zeros(ns+1,n_elements+1)
    for j in 1:ns+1
        eta = (j-1)/ns
        left = min(searchsortedlast(beamNodes,eta),n_elements)
        f = (eta-beamNodes[left])/(beamNodes[left+1]-beamNodes[left])
        weights[j,left],weights[j,left+1] = 1-f,f
    end
    grid,_,_ = wing_geometry(model,chord,sparFraction,nc,weights;sweepAngle=sweep)
    dt = 0.25chord/U
    b = span*cos(sweep)
    reference = UVLM.Reference(b*chord,chord,b,zeros(3),U,1.225)
    aero_model = UVLM.create_UVLMModel(;surfaces=[grid],reference,symmetric=true)
    solver = UVLM.create_UVLMSolver(;coreRadius=1e-3,wakeSheddingFraction=0.1)
    point = UVLM.create_OperatingPoint(;airspeed=U,angleOfAttack=alpha)
    fs = UVLM.freestream(point)
    R = UVLM.body_to_wind(fs)
    Vbody = UVLM.freestream_velocity(fs)
    Vwind = R*Vbody
    isapprox(Vwind,[U,0,0];atol=1e-12) || error("Wind-axis rotation failed")
    problem = UVLM.create_UVLMDynamicProblem(;model=aero_model,aeroSolver=solver,
        operatingPoint=point,timeVector=collect(0:steps).*dt,maximumWakeRows=steps)
    frames = Any[(time=0.0,wake=copy(problem.state.system.wakes[1][1:0,:]),total=0)]
    for step in 1:steps
        UVLM.begin_time_step!(problem)
        UVLM.evaluate_trial!(problem)
        UVLM.commit_time_step!(problem)
        system = problem.state.system
        count = min(visible_rows,system.nwake[1])
        push!(frames,(time=problem.state.timeNow,
            wake=copy(system.wakes[1][1:count,:]),total=system.nwake[1]))
        if step%10 == 0
            println("Accepted step $step; retained $(system.nwake[1]) wake rows")
            flush(stdout)
        end
    end
    serialize(joinpath(output,"computed_frames.jls"),(;grid,frames,alpha,U,dt,nc,ns,R))
    te_root = grid[:,end,1]
    station = ns÷2+1
    te_mid = grid[:,end,station]
    transform(p) = R*(p-te_root)
    transform_mid(p) = R*(p-te_mid)
    positions = [transform(grid[:,i,j]) for i in 1:nc+1,j in 1:ns+1]
    sx,sy,sz = lattice_lines(positions)
    te = [transform(grid[:,end,j]) for j in 1:ns+1]
    tex,tey,tez = ([p[d] for p in te] for d in 1:3)
    function panel_lines(wake)
        x,y,z = Float64[],Float64[],Float64[]
        for panel in wake
            for p in (panel.rtl,panel.rtr,panel.rbr,panel.rbl,panel.rtl)
                q = transform(p)
                push!(x,q[1]);push!(y,q[2]);push!(z,q[3])
            end
            push!(x,NaN);push!(y,NaN);push!(z,NaN)
        end
        return x,y,z
    end
    # Identical metre ranges on x/y/z and aspect_ratio=:equal prevent any
    # vertical exaggeration. Limits remain fixed for the entire animation.
    points = collect(vec(positions))
    for frame in frames,p in frame.wake,corner in (:rtl,:rtr,:rbl,:rbr)
        push!(points,transform(getproperty(p,corner)))
    end
    lows = [minimum(p[d] for p in points) for d in 1:3]
    highs = [maximum(p[d] for p in points) for d in 1:3]
    width = 1.12maximum(highs-lows)
    centers = (lows+highs)/2
    limits = [(centers[d]-width/2,centers[d]+width/2) for d in 1:3]
    side_x = (-0.12,0.30)
    side_z = (-0.21,0.21)
    common = (;legend=false,grid=true,framestyle=:box,
        aspect_ratio=:equal,guidefontsize=10,tickfontsize=8,titlefontsize=11)
    Plots.gr()
    function figure_for(frame)
        wx,wy,wz = panel_lines(frame.wake)
        plan = Plots.plot(sx,sy;common...,color=:steelblue,lw=1.1,
            xlabel="Wind x [m] (freestream direction)",ylabel="Span y [m]",
            xlims=limits[1],ylims=limits[2],title="Plan view | equal metre scale")
        spatial = Plots.plot3d(sx,sy,sz;common...,color=:steelblue,lw=1.0,
            xlabel="Wind x [m]",ylabel="y [m]",zlabel="Wind z [m]",
            xlims=limits[1],ylims=limits[2],zlims=limits[3],camera=(35,25),
            title="3D | equal metre scale on all axes")
        side_grid = [transform_mid(grid[:,i,station]) for i in 1:nc+1]
        side = Plots.plot([p[1] for p in side_grid],[p[3] for p in side_grid];
            common...,color=:steelblue,lw=2.2,xlims=side_x,ylims=side_z,
            xlabel="Wind x aft of midspan TE [m]",ylabel="Wind z relative to TE [m]",
            title="Midspan side view | equal metre scale")
        Plots.hline!(side,[0];color=:gray50,ls=:dash,lw=1)
        if !isempty(frame.wake)
            nx,ny,nz = panel_lines(frame.wake[1:1,:])
            Plots.plot!(plan,wx,wy;color=:darkorange,lw=0.9)
            Plots.plot!(plan,nx,ny;color=:crimson,lw=1.4)
            Plots.plot3d!(spatial,wx,wy,wz;color=:darkorange,lw=0.9)
            Plots.plot3d!(spatial,nx,ny,nz;color=:crimson,lw=1.4)
            for (r,panel) in enumerate(frame.wake[:,station])
                a,b = transform_mid(panel.rtl),transform_mid(panel.rbl)
                Plots.plot!(side,[a[1],b[1]],[a[3],b[3]];
                    color=r==1 ? :crimson : :darkorange,lw=2,marker=:circle,ms=2)
            end
        end
        Plots.plot!(plan,tex,tey;color=:black,lw=2)
        Plots.plot3d!(spatial,tex,tey,tez;color=:black,lw=2)
        Plots.scatter!(side,[0],[0];color=:black,ms=3)
        Plots.annotate!(side,0.17,0.045,Plots.text("U = 45 m/s  →",10))
        Plots.annotate!(side,-0.075,-0.038,Plots.text("wing at 3 deg",9))
        heading = @sprintf("S20 Pazy | AoA = 3 deg | t = %.2f ms | nearest %d of %d wake rows | wind axes; no vertical exaggeration",
            1000frame.time,size(frame.wake,1),frame.total)
        return Plots.plot(plan,spatial,side;layout=(1,3),size=(1800,720),dpi=140,
            plot_title=heading,plot_titlefontsize=14,bottom_margin=5Plots.mm,
            left_margin=13Plots.mm)
    end
    Plots.savefig(figure_for(frames[end]),joinpath(output,"s20_equal_scale.png"))
    Plots.savefig(figure_for(frames[end]),joinpath(output,"s20_equal_scale.pdf"))
    Plots.savefig(figure_for(frames[11]),joinpath(output,"s20_startup_equal_scale.png"))
    animation = Plots.Animation()
    for index in unique(vcat(collect(1:11),collect(16:5:length(frames))))
        Plots.frame(animation,figure_for(frames[index]))
    end
    Plots.gif(animation,joinpath(output,"s20_equal_scale.gif");fps=4)
    loads = UVLM.dimensional_loads(problem)
    force_wind = R*sum(loads.forces[1])
    open(joinpath(output,"verification.txt"),"w") do io
        println(io,"Freestream in body axes [m/s]: $Vbody")
        println(io,"Freestream in wind axes [m/s]: $Vwind")
        println(io,"Final aerodynamic force in wind axes [N]: $force_wind")
        println(io,"3D axis widths [m]: $([b-a for (a,b) in limits])")
        println(io,"Side-view axis widths [m]: $(side_x[2]-side_x[1]), $(side_z[2]-side_z[1])")
        println(io,"Only coordinates transformed; aerodynamic solver unchanged. Fixed undeformed wing.")
        for index in (11,length(frames))
            frame = frames[index]
            println(io,"\nTime=$(frame.time) s; row, x_midspan_mm, z_midspan_mm")
            for (r,panel) in enumerate(frame.wake[:,station])
                p = transform_mid(panel.rbl)
                println(io,"$r,$(1000p[1]),$(1000p[3])")
            end
        end
    end
    println(read(joinpath(output,"verification.txt"),String))
    println("Rendered $output")
end

render_equal_scale_wake()
