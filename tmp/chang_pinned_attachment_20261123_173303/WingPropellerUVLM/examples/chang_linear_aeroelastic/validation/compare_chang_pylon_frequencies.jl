# Compare the two reduced pylon representations with Chang et al. (2024).
# Run from the repository root:
# julia --project=lib/WingPropellerUVLM/studies lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/validation/compare_chang_pylon_frequencies.jl
module ChangPylonFrequencyComparison

using LinearAlgebra,Printf
get!(ENV,"GKSwstype","100")
using Plots

const EXAMPLE_DIR = normpath(joinpath(@__DIR__,".."))
include(joinpath(EXAMPLE_DIR,"src","ChangAeroelastic.jl"))
using .ChangAeroelastic

const OUTPUT_DIR = joinpath(EXAMPLE_DIR,"output","validation","chang_pylon_comparison")
const SOURCE_URL = "https://www.researchgate.net/publication/380374261_Wing_and_Propeller_Aerodynamic_Interaction_Effects_on_Whirl_Flutter_Instability_VFS_80th_Annual_Forum"
const SOURCE_DOI = "https://doi.org/10.4050/F-0080-2024-1330"
# Transcribed from the author-uploaded preprint, Tables 3, 4, and 7.
# Table 7 is explicitly for the flexible wing with the propeller at rest.
const COUPLED_REFERENCE = [
    (mode="Wing OOP1",family=:wing_oop,hz=5.89),
    (mode="Pylon IP1",family=:pylon_ip,hz=7.95),
    (mode="Pylon OOP1",family=:pylon_oop,hz=8.55),
    (mode="Pylon T1",family=:pylon_torsion,hz=12.81),
    (mode="Wing IP1",family=:wing_ip,hz=19.51)
]
const SUPPORT_REFERENCE = [
    (mode="Pylon OOP1",family=:pitch,hz=7.97),
    (mode="Pylon IP1",family=:yaw,hz=7.97),
    (mode="Pylon T1",family=:not_retained,hz=12.73),
    (mode="Pylon OOP2",family=:not_retained,hz=38.17),
    (mode="Pylon IP2",family=:not_retained,hz=38.17)
]
const WING_REFERENCE = [
    (mode="Wing OOP1",family=:out_of_plane,order=1,hz=6.63),
    (mode="Wing IP1",family=:in_plane,order=1,hz=20.93),
    (mode="Wing OOP2",family=:out_of_plane,order=2,hz=36.66),
    (mode="Wing T1",family=:torsion,order=1,hz=43.50),
    (mode="Wing OOP3",family=:out_of_plane,order=3,hz=93.31)
]

function assemble_retuned_rigid(p,config)
    # Same rigid-body closure as output/validation/pylon_comparison_audit.jl:
    # uniform slender pylon; thin axisymmetric rotor; no added nacelle inertia.
    L = p.L_pylon
    mr = p.Nb_prop * config.structural.blade_mass_kg
    mp = config.structural.pylon_mass_per_length_kgpm * L
    S = -(mr * L + mp * L / 2)
    Jt = p.Ix_prop / 2
    Ir = mr * L^2 + mp * L^2 / 3 + Jt
    rigidParameters = merge(
        p,
        (;
            Inθ_prop=Ir,Inψ_prop=Ir,
            Kθ_prop=Ir * (2pi * config.structural.pitch_frequency_hz)^2,
            Kψ_prop=Ir * (2pi * config.structural.yaw_frequency_hz)^2,
            SθP_prop=S,SψP_prop=S,SαP_prop=S,SγP_prop=S,
            IθαP_prop=Ir,IψγP_prop=Ir,IαP_prop=Ir,IγP_prop=Ir
        )
    )
    structural = assemble_structural_model(rigidParameters)
    for node in p.prop_attach_nodes
        fullIndex = (node - 1) * p.ndof + 5
        structural.M_global[fullIndex,fullIndex] += p.Ix_prop
        freeIndex = fullIndex - p.ndof
        freeIndex > 0 && (structural.M[freeIndex,freeIndex] += p.Ix_prop)
    end
    return structural,rigidParameters
end

function coupled_modes(s)
    @assert isposdef(Symmetric(s.M)) && isposdef(Symmetric(s.K))
    @assert isapprox(s.M,s.M'; rtol=1e-12,atol=1e-12)
    @assert isapprox(s.K,s.K'; rtol=1e-12,atol=1e-12)
    @assert norm(s.C,Inf) == 0.0  # Stopped rotor, no damping or aero loads.
    solution = eigen(Symmetric(s.K),Symmetric(s.M))
    n = s.ndof_wing_free
    families = (
        wing_oop=reduce(vcat,(i .+ [3,5] for i = 0:6:(n-6))),
        wing_ip=reduce(vcat,(i .+ [2,6] for i = 0:6:(n-6))),
        wing_torsion=collect(4:6:n),
        wing_axial=collect(1:6:n),
        pylon_oop=[n + 1],
        pylon_ip=[n + 2]
    )
    modes = NamedTuple[]
    maximumResidual = 0.0
    for (rank,value) in enumerate(solution.values)
        @assert isfinite(value) && value > 0
        q = solution.vectors[:,rank]
        # Projected strain-energy fractions identify the dominant family.
        # This follows the existing wing classifier; cross-family strain terms
        # are not allocated. It is not a MAC match to unavailable RCAS vectors.
        energies = map(families) do indices
            max(real(dot(q[indices],s.K[indices,indices] * q[indices])),0.0)
        end
        fractions = map(e -> e / sum(energies),energies)
        family = collect(keys(families))[argmax(collect(values(energies)))]
        residual = norm(s.K * q - value * s.M * q) /
                   ((opnorm(s.K,Inf) + value * opnorm(s.M,Inf)) * norm(q))
        maximumResidual = max(maximumResidual,residual)
        push!(modes,(; rank,frequency_hz=sqrt(value) / (2pi),family,fractions))
    end
    @assert maximumResidual < 1e-10
    return modes,maximumResidual
end

function comparison_row(setting,scope,table,mode,reference,rigid,assumed;
    rigid_rank=missing,assumed_rank=missing)
    error_percent(value) = ismissing(value) ? missing : 100 * (value / reference - 1)
    return (; setting,scope,source_table=table,mode,chang_hz=reference,
        rigid_retuned_hz=rigid,rigid_error_percent=error_percent(rigid),
        assumed_mode_hz=assumed,assumed_error_percent=error_percent(assumed),
        rigid_rank,assumed_rank,
        status=ismissing(rigid) && ismissing(assumed) ? "not_retained" : "represented")
end

csv_cell(x) = ismissing(x) ? "" : "\"" * replace(string(x),"\"" => "\"\"") * "\""
function write_csv(path,rows)
    open(path,"w") do io
        println(io,join(string.(keys(first(rows))),','))
        for row in rows
            println(io,join(csv_cell.(values(row)),','))
        end
    end
end

fmt(x; digits=4) = ismissing(x) ? "Not retained" : @sprintf("%.*f",digits,x)
function write_markdown_table(io,rows)
    println(io,"| Mode | Chang (Hz) | Rigid retuned (Hz) | Error (%) | Assumed mode (Hz) | Error (%) |")
    println(io,"|---|---:|---:|---:|---:|---:|")
    for r in rows
        println(
            io,
            "| $(r.mode) | $(fmt(r.chang_hz;digits=2)) | $(fmt(r.rigid_retuned_hz)) | $(fmt(r.rigid_error_percent;digits=2)) | $(fmt(r.assumed_mode_hz)) | $(fmt(r.assumed_error_percent;digits=2)) |"
        )
    end
end

function plot_comparison(rows,title,stem)
    colors = (:gray30,:darkorange2,:steelblue)
    labels = ["Wing\nOOP1","Pylon\nIP1","Pylon\nOOP1","Pylon\nT1","Wing\nIP1"]
    x = collect(1:length(rows))
    frequencies = plot(; title,ylabel="Natural frequency (Hz)",
        xticks=(x,labels),ylims=(0,23),legend=:topleft,legendfontsize=8,
        grid=:y,foreground_color_legend=nothing,background_color_legend=:transparent)
    columns = (:chang_hz,:rigid_retuned_hz,:assumed_mode_hz)
    legends = ("Chang RCAS (Table 7)","Rigid pylon; retuned stiffness","Assumed-mode pylon")
    for j = 1:3
        y = [ismissing(r[columns[j]]) ? NaN : r[columns[j]] for r in rows]
        bar!(frequencies,x .+ (j - 2) * 0.24,y; bar_width=0.22,
            color=colors[j],linecolor=colors[j],label=legends[j])
    end
    annotate!(frequencies,4.0,15.0,text("Torsion absent\nin both reduced\nmodels",8,:center))
    represented = filter(r -> r.status == "represented",rows)
    errorValues = vcat([r.rigid_error_percent for r in represented],
        [r.assumed_error_percent for r in represented])
    errorLimits = (min(minimum(errorValues),0.0) - 0.5,
        max(maximum(errorValues),0.0) + 0.5)
    errors = plot(; title="Signed error relative to Chang",ylabel="Frequency error (%)",
        xticks=(1:4,labels[[1,2,3,5]]),ylims=errorLimits,legend=false,grid=:y)
    for (j,column) in enumerate((:rigid_error_percent,:assumed_error_percent))
        bar!(errors,(1:4) .+ (j == 1 ? -0.14 : 0.14),
            [r[column] for r in represented]; bar_width=0.26,
            color=colors[j+1],linecolor=colors[j+1])
    end
    hline!(errors,[0.0]; color=:gray30,linewidth=1,label=false)
    figure = plot(frequencies,errors; layout=(1,2),size=(1300,500),
        margin=5 * Plots.mm,top_margin=8 * Plots.mm,bottom_margin=9 * Plots.mm,
        titlefontsize=11,guidefontsize=10,tickfontsize=9,dpi=180)
    savefig(figure,joinpath(OUTPUT_DIR,stem * ".png"))
    savefig(figure,joinpath(OUTPUT_DIR,stem * ".svg"))
end

function main()
    defaults = chang_case_defaults()
    @assert length(defaults.propeller.attachment_eta) == 1 "This comparison requires one propeller."
    aligned = merge(
        defaults,
        (;
            wing=merge(defaults.wing,(; spanwise_panels=30)),
            propeller=merge(defaults.propeller,(; attachment_eta=[6.25 / defaults.wing.span_m]))
        )
    )
    settings = (("current_defaults",defaults),("reference_aligned_30",aligned))
    rows,spectra,metadata = NamedTuple[],NamedTuple[],NamedTuple[]
    for (setting,caseDefaults) in settings
        config = load_chang_configuration(caseDefaults; env=Dict{String,String}())
        p = ChangAeroelastic.build_chang_model_parameters(config)
        p = merge(p,(; Ω=0.0,ξ_prop=0.0,stiffness_damping_ratio=0.0))
        assumed = assemble_structural_model(p)
        rigid,rp = assemble_retuned_rigid(p,config)
        am,ar = coupled_modes(assumed)
        rm,rr = coupled_modes(rigid)
        station = p.span_nodes[only(p.prop_attach_nodes)]
        setting == "reference_aligned_30" && @assert isapprox(station,6.25; atol=1e-12)
        push!(
            metadata,
            (; setting,elements=p.Ne,attachment_node=only(p.prop_attach_nodes),
                requested_station_m=only(p.propeller_span_positions),actual_station_m=station,
                assumed_pitch_stiffness=p.Kθ_prop,assumed_yaw_stiffness=p.Kψ_prop,
                rigid_pitch_stiffness=rp.Kθ_prop,rigid_yaw_stiffness=rp.Kψ_prop,
                rigid_pitch_inertia=rp.Inθ_prop,assumed_pitch_inertia=p.Inθ_prop,
                assumed_yaw_inertia=p.Inψ_prop,max_eigen_residual=max(ar,rr))
        )
        for (model,modes) in (("rigid_retuned",rm),("assumed_mode",am))
            for m in modes
                push!(
                    spectra,
                    merge((; setting,model,rank=m.rank,
                            frequency_hz=m.frequency_hz,dominant_family=m.family),m.fractions)
                )
            end
        end
        for ref in COUPLED_REFERENCE
            if ref.family == :pylon_torsion
                push!(rows,comparison_row(setting,"coupled",7,ref.mode,ref.hz,missing,missing))
                continue
            end
            ai = findfirst(m -> m.family == ref.family,am)
            ri = findfirst(m -> m.family == ref.family,rm)
            isnothing(ai) && error("No assumed-mode match for $(ref.mode)")
            isnothing(ri) && error("No rigid-model match for $(ref.mode)")
            push!(
                rows,
                comparison_row(setting,"coupled",7,ref.mode,ref.hz,
                    rm[ri].frequency_hz,am[ai].frequency_hz;
                    rigid_rank=rm[ri].rank,assumed_rank=am[ai].rank)
            )
        end
        # Both representations share the same bare-wing model. Do not compare
        # these isolated frequencies to the installed-system modes of Table 7.
        wing_modes = ChangAeroelastic.chang_wing_modal_analysis(assumed; number_of_modes=10)
        for ref in WING_REFERENCE
            index = findfirst(m -> m.family == ref.family && m.family_order == ref.order,wing_modes)
            isnothing(index) && error("No wing-only match for $(ref.mode)")
            hz = wing_modes[index].frequency_hz
            push!(rows,comparison_row(setting,"isolated_wing",4,ref.mode,ref.hz,hz,hz))
        end
        if setting == "current_defaults"
            for ref in SUPPORT_REFERENCE
                a =
                    ref.family == :pitch ? sqrt(p.Kθ_prop / p.Inθ_prop) / (2pi) :
                    ref.family == :yaw ? sqrt(p.Kψ_prop / p.Inψ_prop) / (2pi) : missing
                r =
                    ref.family == :pitch ? sqrt(rp.Kθ_prop / rp.Inθ_prop) / (2pi) :
                    ref.family == :yaw ? sqrt(rp.Kψ_prop / rp.Inψ_prop) / (2pi) : missing
                push!(rows,comparison_row(setting,"isolated_support",3,ref.mode,ref.hz,r,a))
            end
        end
    end
    mkpath(OUTPUT_DIR)
    write_csv(joinpath(OUTPUT_DIR,"frequency_comparison.csv"),rows)
    write_csv(joinpath(OUTPUT_DIR,"full_nonspinning_spectra.csv"),spectra)
    write_csv(joinpath(OUTPUT_DIR,"model_settings.csv"),metadata)
    report_path = joinpath(OUTPUT_DIR,"comparison.md")
    open(report_path,"w") do io
        println(io,"# Chang RCAS, retuned rigid pylon, and assumed-mode pylon\n")
        println(
            io,
            "Reference: Chang, Sanghi, and Cesnik, *Wing and Propeller Aerodynamic Interaction Effects on Whirl Flutter Instability*, VFS Forum 80 (2024). [DOI]($SOURCE_DOI); [author preprint]($SOURCE_URL). Tables 3, 4, and 7 were transcribed from the preprint.\n"
        )
        println(
            io,
            "All comparisons use a stopped rotor, no aerodynamic loading, and no structural damping. Errors are 100*(computed/Chang - 1). Table 7 concerns the installed flexible-wing system; Tables 3 and 4 concern isolated support and isolated wing respectively.\n"
        )
        println(
            io,
            "The rigid reference uses a full-length uniform slender pylon, a thin axisymmetric rotor with transverse inertia Ix/2, no additional nacelle mass, and rotor axial inertia carried by the wing. Its stiffnesses are retuned to the configured isolated pitch/yaw targets; they are not Chang's adjusted-pylon stiffnesses. The assumed-mode production formulation is unchanged.\n"
        )
        println(
            io,
            "Coupled modes are matched by their dominant projected strain-energy family, using pitch for pylon OOP and yaw for pylon IP. Cross-family energy terms are not allocated. RCAS mode vectors are unavailable, so this is a family correspondence, not a MAC validation. Missing torsion is left blank in CSVs and omitted from error calculations.\n"
        )
        for data in metadata
            println(io,"## $(data.setting)\n")
            println(
                io,
                "Wing elements: $(data.elements). Requested attachment: $(data.requested_station_m) m; actual attachment: $(data.actual_station_m) m at node $(data.attachment_node). Chang specifies 6.25 m. The aligned case changes only mesh and attachment setting; it retains the production pylon length and other inputs.\n"
            )
            println(
                io,
                "Pitch/yaw stiffnesses (N m/rad): assumed $(data.assumed_pitch_stiffness) / $(data.assumed_yaw_stiffness); rigid retuned $(data.rigid_pitch_stiffness) / $(data.rigid_yaw_stiffness).\n"
            )
            println(io,"Installed system, Chang Table 7:\n")
            selected = filter(r -> r.setting == data.setting && r.scope == "coupled",rows)
            write_markdown_table(io,selected)
            println(io,"\nBare wing, Chang Table 4 (identical wing model in both columns):\n")
            write_markdown_table(io,filter(r -> r.setting == data.setting && r.scope == "isolated_wing",rows))
            println(
                io,
                "\nMass and stiffness are symmetric positive definite. Maximum relative eigenproblem residual: $(data.max_eigen_residual).\n"
            )
            title =
                data.setting == "current_defaults" ?
                "Current defaults: $(data.elements) elements; attachment $(data.actual_station_m) m" :
                "Reference station: 30 elements; attachment 6.25 m"
            plot_comparison(selected,title,data.setting)
            println(io,"![Frequency and error comparison]($(data.setting).png)\n")
        end
        println(io,"## Isolated support, Chang Table 3\n")
        write_markdown_table(io,filter(r -> r.scope == "isolated_support",rows))
        println(
            io,
            "\nThe 7.97 Hz values are calibration targets in both reduced models. Neither model retains pylon torsion or second pylon bending modes. The twist frequency/stiffness pair in the production inputs defines rotor axial inertia; it does not add a torsional coordinate.\n"
        )
        println(
            io,
            "All computed spectra and dominant-family fractions are saved in full_nonspinning_spectra.csv. Settings are saved in model_settings.csv. Reproduce from the repository root:\n\n```powershell\njulia --project=lib/WingPropellerUVLM/studies lib/WingPropellerUVLM/examples/chang_linear_aeroelastic/validation/compare_chang_pylon_frequencies.jl\n```\n"
        )
    end
    for data in metadata
        println("\n",data)
        write_markdown_table(stdout,filter(r -> r.setting == data.setting && r.scope == "coupled",rows))
    end
    println("\nSaved comparison: ",report_path)
    return (; rows,spectra,metadata)
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    ChangPylonFrequencyComparison.main()
end
