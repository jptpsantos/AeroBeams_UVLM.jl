# Run directly with julia --project=lib/WingPropellerUVLM/studies \
#     lib/WingPropellerUVLM/test/chang_aeroelastic_convergence_plots.jl
module ChangAeroelasticConvergencePlotTests
using Test,DelimitedFiles,TOML

const EXAMPLE = normpath(joinpath(@__DIR__,"..","examples","chang_linear_aeroelastic"))
include(joinpath(EXAMPLE,"studies","convergence","plot_chang_convergence_aeroelastic.jl"))
const P = ChangConvergenceAeroelasticPlots
const ARCHIVES = joinpath(EXAMPLE,"output","convergence_aeroelastic_explicit")
const CASE_NAMES = ("First case","Second case")

function copy_tables(directory)
    for name in CASE_NAMES
        target = joinpath(directory,name)
        mkpath(target)
        for file in ("study.toml","case_matrix.csv","convergence.csv")
            cp(joinpath(ARCHIVES,name,file),joinpath(target,file))
        end
    end
    return directory
end

function edit_result!(edit!,directory; case="Second case",family="prop_radial",level=1)
    path = joinpath(directory,case,"convergence.csv")
    data,header = readdlm(path,','; header=true)
    columns = Dict(string(name) => i for (i,name) in enumerate(vec(header)))
    row = only(
        i for i in axes(data,1) if string(data[i,columns["family"]]) == family &&
        parse(Float64,string(data[i,columns["level"]])) == level
    )
    edit!(data,row,columns)
    open(path,"w") do io
        writedlm(io,header,',')
        writedlm(io,data,',')
    end
end

@testset "Combined archived Chang blade damping convergence" begin
    result = P.load_results(ARCHIVES)
    @test result.input_count == 16
    @test length(result.cases) == 11
    @test length(result.rows) == 12
    @test Set(r.family for r in result.rows) == Set(("prop_radial","prop_chord"))
    @test length(unique(c.id for c in result.cases)) == 11
    @test count(c -> length(c.families) == 2,result.cases) == 1

    radial = sort(filter(r -> r.family == "prop_radial",result.rows); by=r -> r.x)
    chord = sort(filter(r -> r.family == "prop_chord",result.rows); by=r -> r.x)
    @test radial[2].case_id == chord[2].case_id
    @test radial[2].values[1:2] == (-1.247861507734489,-1.2488302066799508)
    for group in (radial,chord)
        @test [r.x for r in group] == [5.0,10.0,20.0,30.0,40.0,50.0]
        @test all(isnan,first(group).previous)
        @test all(iszero,last(group).errors)
        @test all(r -> r.reference_panels == 50.0,group)
        @test all(r -> r.limits == (0.01,0.01,0.1,0.1),group)
    end

    # The archives use different references (20 and 50 panels). These fixed
    # values catch accidentally concatenating their existing error columns.
    @test radial[1].errors[1] ≈ 0.130594687107213
    @test radial[3].errors[1] ≈ 0.0186369851137851
    @test radial[4].previous[1] ≈ 0.0127190924276654
    @test radial[5].previous[2] ≈ 0.00366718863159599
    @test chord[1].errors[2] ≈ 0.0612680082776527
    @test chord[3].errors[1] ≈ 0.0112702172787922
    @test chord[4].previous[2] ≈ 0.00537583504681094
    @test radial[end].values[1:2] == (-1.1893704311219957,-1.1902486247094397)
    @test chord[end].values[1:2] == (-1.2724227945957598,-1.2733675995637486)

    # Only small saved tables are needed; histories and simulations are unused.
    mktempdir() do directory
        copy_tables(directory)
        copied = P.load_results(directory)
        @test [r.values for r in copied.rows] == [r.values for r in result.rows]
    end
end

@testset "Damping errors in percent" begin
    result = P.load_results(ARCHIVES)
    radial30 = only(r for r in result.rows if r.family == "prop_radial" && r.x == 30)
    radial40 = only(r for r in result.rows if r.family == "prop_radial" && r.x == 40)
    chord30 = only(r for r in result.rows if r.family == "prop_chord" && r.x == 30)

    # Fixed numerical values catch a missing factor of 100, a signed denominator,
    # or use of the previous/finest damping rate as the successive denominator.
    @test isapprox(P.successive_errors_percent(radial30)[1],1.064102457484329)
    @test isapprox(P.successive_errors_percent(radial40)[2],0.3075140861314764)
    @test isapprox(P.successive_errors_percent(chord30)[2],0.4241427465050588)
    @test isapprox(P.relative_limits_percent(radial30)[1],0.8366182284907304)
    @test P.relative_error_percent(0.1,-2.0) == 5.0
    @test P.relative_error_percent(0.0,-2.0) == 0.0
    @test isnan(P.relative_error_percent(0.1,0.0))
    @test all(isnan,P.successive_errors_percent(first(result.rows)))
    radial = filter(r -> r.family == "prop_radial",result.rows)
    @test P.reference_errors_percent(first(radial),last(radial))[1] ≈ 10.980152498327723
    @test P.reference_errors_percent(radial30,last(radial))[1] ≈ 0.4975651429754296
    @test P.reference_errors_percent(last(radial),last(radial)) == (0.0,0.0)
    zero_reference = P.reference_errors_percent((; errors=(0.1,0.2)),(; values=(0.0,-2.0)))
    @test isnan(zero_reference[1])
    @test zero_reference[2] == 10.0

    mktempdir() do directory
        P.write_tables(directory,result)
        saved = P.read_table(joinpath(directory,"damping_convergence.csv"))
        row(family,panels) = only(r for r in saved if r["family"] == family && r["panels"] == panels)
        @test isapprox(row("prop_radial",30)["previous_relative_pitch_percent"],1.064102457484329)
        @test isapprox(row("prop_radial",40)["previous_relative_yaw_percent"],0.3075140861314764)
        @test isapprox(row("prop_chord",30)["previous_relative_yaw_percent"],0.4241427465050588)
        @test isnan(row("prop_radial",5)["previous_relative_pitch_percent"])
        @test row("prop_chord",50)["reference_relative_yaw_percent"] == 0.0
        @test row("prop_radial",5)["reference_relative_pitch_percent"] ≈ 10.980152498327723
    end
end

@testset "Reject incompatible saved damping results" begin
    mktempdir() do directory
        copy_tables(directory)
        path = joinpath(directory,"Second case","study.toml")
        study = TOML.parsefile(path)
        study["settings"]["physical"]["alpha_deg"] = 4.0
        open(io -> TOML.print(io,study),path,"w")
        @test_throws ErrorException P.load_results(directory)
    end

    mktempdir() do directory
        copy_tables(directory)
        # This is the same 10-by-10 baseline as in the first archive.
        edit_result!(directory) do data,row,columns
            data[row,columns["pitch_lambda_per_s"]] = -1.0
        end
        @test_throws ErrorException P.load_results(directory)
    end

    mktempdir() do directory
        copy_tables(directory)
        edit_result!(directory; level=5) do data,row,columns
            data[row,columns["status"]] = "failed"
        end
        @test_throws ErrorException P.load_results(directory)
    end

    mktempdir() do directory
        copy_tables(directory)
        edit_result!(directory; family="prop_chord",level=5) do data,row,columns
            data[row,columns["valid"]] = false
        end
        @test_throws ErrorException P.load_results(directory)
    end
end
end
