using Test,DelimitedFiles

const EXAMPLE = normpath(joinpath(@__DIR__,".."))

@testset "Primary runner and dependency-free organized exports" begin
    mktempdir() do directory
        withenv(
            "CHANG_WING_SPAN_PANELS" => "4","CHANG_WING_CHORD_PANELS" => "2",
            "CHANG_PROP_RADIAL_PANELS" => "2","CHANG_PROP_CHORD_PANELS" => "2",
            "CHANG_WAKE_ROWS_WING" => "2","CHANG_WAKE_ROWS_PROPELLER" => "2",
            "CHANG_END_TIME_S" => "0.001","CHANG_PLOT_RESULTS" => "false",
            "CHANG_ANIMATE_WAKE" => "false","CHANG_OUTPUT_DIR" => directory
        ) do
            include(joinpath(EXAMPLE,"run_chang_linear_aeroelastic.jl"))
        end
        @test chang_run.results.validation_passed
        single = Base.invokelatest(export_chang_organized_history,chang_run)
        @test dirname(single.csv_path) == joinpath(directory,"multi_prop")
        data,headers = readdlm(single.csv_path,',',Float64; header=true)
        @test vec(headers) == string.(collect(keys(single.df)))
        @test data[:,1] == chang_run.results.time
        @test data[:,4] == rad2deg.(chang_run.results.propeller_pitch[1])
        @test isfile(single.txt_path)

        # Distinct synthetic histories make the inboard/outboard ordering observable.
        config = merge(chang_run.config,(; propeller=merge(chang_run.config.propeller,
            (; attachment_eta=[0.8,0.4]))))
        n = length(chang_run.results.time)
        results = merge(
            chang_run.results,
            (;
                propeller_pitch=[fill(0.1,n),fill(0.2,n)],
                propeller_yaw=[fill(0.3,n),fill(0.4,n)])
        )
        twin = Base.invokelatest(export_chang_organized_history,(; config,results))
        data,headers = readdlm(twin.csv_path,',',Float64; header=true)
        @test vec(headers) == string.(collect(keys(twin.df)))
        @test data[:,4] == fill(rad2deg(0.2),n)
        @test data[:,6] == fill(rad2deg(0.1),n)
        @test occursin("inboard_propeller_index = 2",read(twin.txt_path,String))
        @test occursin("outboard_propeller_index = 1",read(twin.txt_path,String))
        @test isfile(single.csv_path)
    end
end
