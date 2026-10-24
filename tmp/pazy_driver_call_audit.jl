# Execute original/current orchestration with recording solver stubs. This
# compares all solver inputs at production ramp durations without solving.
using Test
const BASE = normpath(joinpath(@__DIR__,"..","lib","WingPropellerUVLM","examples","pazy_aerobeams_coupled"))
function load_driver(name,path)
    scope = Module(name)
    Core.eval(scope,quote
        using DelimitedFiles
        const recorded_static = Any[]
        const recorded_dynamic = Any[]
        const Plots = (;plot=(args...;kwargs...)->nothing,
            plot! = (args...;kwargs...)->nothing,twinx=x->nothing,
            savefig=(args...)->nothing,mm=1.)
        const UVLM = (;create_UVLMSolver=(;kwargs...)->nothing)
        const AeroBeams = (;geometrical_properties_Pazy=()->(15,.549843728,.0989,.3),
            gcd_floats_approx=x->Float64(gcd(Int.(x)...)))
        pazy_lco_branch_names(s;option_name="plot_branches") =
            s == :both ? ("upsweep","downsweep") : s == :upsweep ? ("upsweep",) : ("downsweep",)
        pazy_LCO_strain_plot(args...;kwargs...) = nothing
        function run_pazy_wing_uvlm_static(;kwargs...)
            push!(recorded_static,(;kwargs...))
            return (;history_columns=["dummy"],history=zeros(1,1),
                converged=true,termination_reason="recording stub",geometry_residual=0.,load_residual=0.)
        end
        function run_pazy_wing_uvlm(;kwargs...)
            push!(recorded_dynamic,(;kwargs...))
            time=[0.,kwargs[:duration]/2,kwargs[:duration]]
            return (;time,airspeed_history=kwargs[:airspeed_profile].(time),
                root_axial_strain=zeros(3),root_curvature_2=zeros(3),root_curvature_3=zeros(3),
                tip_twist_degrees=zeros(3),tip_out_of_plane=zeros(3))
        end
    end)
    source = read(path,String)
    # Remove only dependency-loading preamble. Execute the actual functions.
    start = first(findfirst("\"\"\"One linear ramp",source))
    Base.include_string(scope,source[start:end],path)
    return scope
end
old = load_driver(:OriginalDriver,joinpath(@__DIR__,"pazy_original_driver_audit.jl"))
new = load_driver(:CurrentDriver,joinpath(BASE,"PazyLinearLCOSweep.jl"))
mktempdir() do dir
    Base.invokelatest(old.run_pazy_linear_LCO;output_directory=joinpath(dir,"old"))
    Base.invokelatest(new.run_pazy_linear_LCO;downsweep_pulse_amplitude=.1,
        downsweep_pulse_start=3.5294117647058827,downsweep_pulse_duration=.1,
        output_directory=joinpath(dir,"new"))
    @test old.recorded_static[1] == new.recorded_static[1]
    a,b = old.recorded_dynamic[1],new.recorded_dynamic[1]
    @test keys(a) == keys(b)
    for key in keys(a)
        key == :airspeed_profile && continue
        @test a[key] == b[key]
    end
    t=range(0,40;length=10001)
    @test Base.invokelatest.(Ref(a.airspeed_profile),t) == Base.invokelatest.(Ref(b.airspeed_profile),t)
    @test a.perturbation_amplitude == b.perturbation_amplitude == 0.
    @test new.recorded_dynamic[2].perturbation_amplitude == .1
    println("PASS: original HEAD and current pulsed driver have identical upsweep static inputs and all ",length(keys(a))," dynamic keyword values (profile checked at 10001 times).")
    println("PASS: current downsweep alone receives pulse coefficient .1 N.")
end
