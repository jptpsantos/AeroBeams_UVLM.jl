ENV["GKSwstype"] = "100"
include(joinpath(@__DIR__, "..", "lib", "WingPropellerUVLM", "examples",
    "pazy_aerobeams_coupled", "PazyWingUVLMCoupling.jl"))
include(joinpath(@__DIR__, "..", "lib", "WingPropellerUVLM", "examples",
    "pazy_aerobeams_coupled", "PazyWingUVLMVisualization.jl"))
using Printf

function render_s20_wake()
    output = joinpath(@__DIR__, "s20_wake_45ms_3deg")
    mkpath(output)
    sweep = deg2rad(20.0)
    alpha = deg2rad(3.0)
    U, nc, ns, nsteps = 45.0, 10, 20, 10
    pazy = create_pazy_uvlm_structure(; sweepAngle=sweep)
    (; model, n_elements, span, chord, sparFraction, beamNodes) = pazy
    weights = zeros(ns + 1, n_elements + 1)
    for j in 1:ns+1
        eta = (j - 1) / ns
        left = min(searchsortedlast(beamNodes, eta), n_elements)
        f = (eta - beamNodes[left]) / (beamNodes[left+1] - beamNodes[left])
        weights[j,left], weights[j,left+1] = 1-f, f
    end
    grid, _, _ = wing_geometry(model, chord, sparFraction, nc, weights;
        sweepAngle=sweep)
    dt = 0.25 * chord / U
    projected_span = span * cos(sweep)
    reference = UVLM.Reference(projected_span * chord, chord, projected_span,
        zeros(3), U, 1.225)
    aero_model = UVLM.create_UVLMModel(; surfaces=[grid], reference,
        symmetric=true, name="S20 wake visualization")
    solver = UVLM.create_UVLMSolver(; coreRadius=1e-3,
        wakeSheddingFraction=0.1, interaction=true)
    point = UVLM.create_OperatingPoint(; airspeed=U, angleOfAttack=alpha)
    problem = UVLM.create_UVLMDynamicProblem(; model=aero_model,
        aeroSolver=solver, operatingPoint=point,
        timeVector=collect(0:nsteps) .* dt, maximumWakeRows=nsteps)
    frames = Any[]
    push!(frames, (time=0.0, wake=copy(problem.state.system.wakes[1][1:0,:])))
    for step in 1:nsteps
        UVLM.begin_time_step!(problem)
        UVLM.evaluate_trial!(problem)
        UVLM.commit_time_step!(problem)
        system = problem.state.system
        wake = copy(system.wakes[1][1:system.nwake[1],:])
        push!(frames, (time=problem.state.timeNow, wake=wake))
        println("Accepted step $step: $(size(wake)) wake panels")
        flush(stdout)
    end

    # Plot stored quadrilateral vertices, preserving the actual shedding gap.
    function panel_lines(wake)
        x,y,z = Float64[],Float64[],Float64[]
        for panel in wake
            for p in (panel.rtl,panel.rtr,panel.rbr,panel.rbl,panel.rtl)
                push!(x,p[1]); push!(y,p[2]); push!(z,p[3])
            end
            push!(x,NaN); push!(y,NaN); push!(z,NaN)
        end
        return x,y,z
    end
    positions = [UVLM.SVector{3}(grid[:,i,j]) for i in 1:nc+1,j in 1:ns+1]
    sx,sy,sz = lattice_lines(positions)
    tex,tey,tez = grid[1,end,:],grid[2,end,:],grid[3,end,:]
    wake = frames[end].wake
    wx,wy,wz = panel_lines(wake)
    nx,ny,nz = panel_lines(wake[1:1,:])
    Plots.gr()
    common = (; legend=false, grid=true, framestyle=:box,
        guidefontsize=10, tickfontsize=8, titlefontsize=11)
    plan = Plots.plot(sx,sy; common..., color=:steelblue, lw=0.9,
        xlabel="Downstream x [m]", ylabel="Span y [m]", aspect_ratio=:equal,
        title="Plan view: 20 spanwise divisions, 10 wake rows")
    Plots.plot!(plan,wx,wy;color=:darkorange,lw=0.8)
    Plots.plot!(plan,nx,ny;color=:crimson,lw=1.4)
    Plots.plot!(plan,tex,tey;color=:black,lw=2)
    spatial = Plots.plot3d(sx,sy,sz;common...,color=:steelblue,lw=0.8,
        xlabel="x [m]",ylabel="y [m]",zlabel="z [m]",camera=(40,25),
        title="Actual UVLM vertices; undeformed wing held fixed")
    Plots.plot3d!(spatial,wx,wy,wz;color=:darkorange,lw=0.8)
    Plots.plot3d!(spatial,nx,ny,nz;color=:crimson,lw=1.4)
    Plots.plot3d!(spatial,tex,tey,tez;color=:black,lw=2)

    station = ns ÷ 2 + 1
    side = Plots.plot(grid[1,:,station],grid[3,:,station]; common...,
        color=:steelblue,lw=2,xlabel="Downstream x [m]",ylabel="z [m]",
        title="Side view at midspan (vertical scale expanded)")
    for (r,panel) in enumerate(wake[:,station])
        Plots.plot!(side,[panel.rtl[1],panel.rbl[1]],
            [panel.rtl[3],panel.rbl[3]];color=r==1 ? :crimson : :darkorange,
            lw=1.7,marker=:circle,ms=2)
    end
    te = grid[:,end,station]
    Plots.scatter!(side,[te[1]],[te[3]];color=:black,ms=4)
    xend = maximum(filter(isfinite,wx))
    Plots.plot!(side,[te[1],xend],[te[3],te[3]+(xend-te[1])*tan(alpha)];
        color=:gray50,ls=:dash,lw=1.1)

    zoom = Plots.plot(;common...,xlabel="Distance aft of TE [mm]",
        ylabel="Height relative to TE [mm]",title="New row at TE: shedding offset and convection",
        xlims=(-2,1000*(wake[3,station].rbl[1]-te[1])+14))
    for r in 1:3
        p = wake[r,station]
        Plots.plot!(zoom,1000 .* ([p.rtl[1],p.rbl[1]] .- te[1]),
            1000 .* ([p.rtl[3],p.rbl[3]] .- te[3]);
            color=r==1 ? :crimson : :darkorange,lw=2,marker=:circle,ms=3)
        Plots.annotate!(zoom,1000*(p.rbl[1]-te[1]),1000*(p.rbl[3]-te[3]),
            Plots.text("  row $r",9,:left))
    end
    Plots.scatter!(zoom,[0.0],[0.0];color=:black,ms=4)
    heading = @sprintf("S20 Pazy | U = 45 m/s | AoA = 3 deg | t = %.3f ms | dt = %.3f ms",1000*nsteps*dt,1000*dt)
    figure = Plots.plot(plan,spatial,side,zoom;layout=(2,2),size=(1600,1050),
        dpi=150,plot_title=heading,plot_titlefontsize=16,
        bottom_margin=5Plots.mm,left_margin=5Plots.mm)
    Plots.savefig(figure,joinpath(output,"s20_wake_actual.png"))
    Plots.savefig(figure,joinpath(output,"s20_wake_actual.pdf"))

    animation = Plots.Animation()
    xlimits = (minimum(grid[1,:,:])-0.015,maximum(filter(isfinite,wx))+0.015)
    ylimits = (-0.01,maximum(grid[2,:,:])+0.01)
    for frame in frames
        p = Plots.plot(sx,sy;common...,color=:steelblue,lw=1,
            xlabel="Downstream x [m]",ylabel="Span y [m]",aspect_ratio=:equal,
            xlims=xlimits,ylims=ylimits,size=(1000,800),
            title=@sprintf("S20 | 45 m/s | AoA 3 deg | t = %.3f ms | %d wake rows",1000*frame.time,size(frame.wake,1)))
        if !isempty(frame.wake)
            ax,ay,_ = panel_lines(frame.wake)
            bx,by,_ = panel_lines(frame.wake[1:1,:])
            Plots.plot!(p,ax,ay;color=:darkorange,lw=0.8)
            Plots.plot!(p,bx,by;color=:crimson,lw=1.4)
        end
        Plots.plot!(p,tex,tey;color=:black,lw=2)
        q = Plots.plot(grid[1,:,station],grid[3,:,station];common...,
            color=:steelblue,lw=2,xlims=xlimits,ylims=(-0.001,0.018),
            xlabel="Downstream x [m]",ylabel="Height z [m]",
            title="Midspan side view (vertical scale expanded)")
        for (r,panel) in enumerate(frame.wake[:,station])
            Plots.plot!(q,[panel.rtl[1],panel.rbl[1]],
                [panel.rtl[3],panel.rbl[3]];
                color=r==1 ? :crimson : :darkorange,lw=2,marker=:circle,ms=2)
        end
        Plots.plot!(q,[te[1],xend],[te[3],te[3]+(xend-te[1])*tan(alpha)];
            color=:gray50,ls=:dash,lw=1.1)
        Plots.scatter!(q,[te[1]],[te[3]];color=:black,ms=3)
        paired = Plots.plot(p,q;layout=(1,2),size=(1500,750),
            bottom_margin=5Plots.mm,left_margin=5Plots.mm)
        Plots.frame(animation,paired)
    end
    Plots.gif(animation,joinpath(output,"s20_wake_shedding.gif");fps=3)
    open(joinpath(output,"wake_vertices.csv"),"w") do io
        println(io,"row,span_panel,corner,x_m,y_m,z_m")
        for r in axes(wake,1),j in axes(wake,2),corner in (:rtl,:rtr,:rbl,:rbr)
            p = getproperty(wake[r,j],corner)
            println(io,join((r,j,corner,p...),','))
        end
    end
    fs = UVLM.freestream_velocity(UVLM.freestream(point))
    open(joinpath(output,"parameters.txt"),"w") do io
        println(io,"Wing held fixed in its undeformed native AeroBeams geometry.")
        println(io,"sweep=20deg alpha=3deg U=45m/s symmetric=true nc=10 ns=20 steps=10")
        println(io,"chord=$chord span=$span dt=$dt eta=0.1 core=0.001")
        println(io,"freestream_velocity=$fs nominal_row_length=$(U*dt)")
        println(io,"shedding_offset=$(0.1 .* fs .* dt)")
        println(io,"wake_z_extrema=$(extrema(filter(isfinite,wz)))")
        println(io,"structural_main_wing_elements=$n_elements")
    end
    println(read(joinpath(output,"parameters.txt"),String))
    println("Rendered $output")
end

render_s20_wake()
