# Read-only model/history audit; no simulation outputs or case settings are changed.
# Run with --project=lib/WingPropellerUVLM. Optional argument: history CSV path.
using LinearAlgebra, DelimitedFiles, Statistics, StaticArrays
using WingPropellerUVLM: imperial_nodal_forces, imperial_nodal_positions
const EXAMPLE = normpath(joinpath(@__DIR__, ".."))
include(joinpath(EXAMPLE, "src", "ChangAeroelastic.jl"))
using .ChangAeroelastic
BLAS.set_num_threads(1)

function main()
    config = load_chang_configuration(; env=Dict{String,String}())
    model = build_chang_model(config)
    p, s = model.parameters, model.structural
    p.Npropellers == 2 || error("This audit requires exactly two propellers")
    println("AUDIT requested_eta=", config.propeller.attachment_eta,
        " actual_eta=", p.span_nodes[p.prop_attach_nodes] ./ p.span_length)
    println("AUDIT dt=", p.dt[1], " wake_rows=", (config.wake.maximum_rows_wing, config.wake.maximum_rows_propeller),
        " nominal_wake_lengths=", p.Vinf * p.dt[1] .* [config.wake.maximum_rows_wing, config.wake.maximum_rows_propeller])
    println("AUDIT symmetry M,K,gyro=", (norm(s.M-s.M'), norm(s.K-s.K'), norm(s.C+s.C')),
        " min_mass_eigenvalue=", eigmin(Symmetric(s.M)))
    reversed_config = merge(config, (; propeller=merge(config.propeller, (;attachment_eta=reverse(config.propeller.attachment_eta)))))
    reversed_model = build_chang_model(reversed_config)
    nw = s.ndof_wing_free
    perm = vcat(1:nw, nw .+ [3,4,1,2])
    println("AUDIT propeller_permutation_errors=", [norm(getproperty(s,key)[perm,perm] - getproperty(reversed_model.structural,key)) for key in (:M,:C,:K)])
    n = size(s.M,1)
    eig = eigen([zeros(n,n) Matrix{Float64}(I,n,n); -(s.M\s.K) -(s.M\s.C)])
    indices = sort(filter(i -> imag(eig.values[i]) > 1e-6, eachindex(eig.values)); by=i -> imag(eig.values[i]))
    for i in indices[1:7]
        v = eig.vectors[1:n,i]
        println("AUDIT gyro_mode frequency_hz=", imag(eig.values[i])/(2pi),
            " sigma=",real(eig.values[i]), " outboard_inboard_angle_ratio=", norm(v[nw+3:nw+4])/norm(v[nw+1:nw+2]))
    end
    # Frozen synthetic dimensional forces isolate load-transfer/indexing errors
    # from aerodynamic discretization and wake convergence.
    w = build_chang_workspace(model)
    w.system.near_field_analysis[] = true
    for arrays in (w.system.span_seg_forces, w.system.chord_seg_forces, w.system.unsteady_forces)
        for (isurf, a) in enumerate(arrays), i in eachindex(a)
            a[i] = SVector(sin(i+isurf), cos(2i-isurf), sin(3i+2isurf))
        end
    end
    forces = imperial_nodal_forces(w.system)
    for amplitude in (0.0, 0.01)
        q = amplitude .* sin.(collect(1:n))
        kin = update_aero_geometry_for_state!(model,w,q,0.02)
        mapped = assemble_structural_aero_load!(model,w,kin)
        errors = Float64[]
        for j in 1:n
            qp, qm = copy(q), copy(q)
            qp[j] += 1e-6
            qm[j] -= 1e-6
            update_aero_geometry_for_state!(model,w,qp,0.02)
            xp = imperial_nodal_positions(w.system)
            update_aero_geometry_for_state!(model,w,qm,0.02)
            xm = imperial_nodal_positions(w.system)
            fd = sum(dot(forces[k][i], (xp[k][i]-xm[k][i])/2e-6) for k in eachindex(forces) for i in eachindex(forces[k]))
            push!(errors, abs(mapped[j]-fd)/max(1.0,abs(mapped[j]),abs(fd)))
        end
        println("AUDIT virtual_work amplitude=",amplitude," maximum_error=",maximum(errors)," prop_errors=",errors[nw+1:end])
    end
    path = isempty(ARGS) ? joinpath(EXAMPLE,"output","chang_linear_imperial_uvlm_history.csv") : abspath(ARGS[1])
    raw, header = readdlm(path, ',', header=true)
    headers = vec(header)
    time = Float64.(raw[:,1])
    for name in ("propeller_1_pitch_deg", "propeller_1_yaw_deg", "propeller_2_pitch_deg", "propeller_2_yaw_deg")
        y = Float64.(raw[:,findfirst(==(name),headers)])
        println("HISTORY ",name," rms_1s_windows=",[sqrt(mean(abs2,y[findall(t-> a<=t<a+1,time)])) for a in 1:4])
        for start in (1.0,2.0,3.0)
            ix = findall(>=(start),time)
            amp = abs.(y .- mean(y[ix]))
            peaks = [j for j in first(ix)+1:last(ix)-1 if amp[j]>amp[j-1] && amp[j]>=amp[j+1] && amp[j]>1e-10]
            x = time[peaks]; z = log.(amp[peaks]); X = hcat(ones(length(x)),x)
            coeff = X\z
            r2 = 1-sum(abs2,z-X*coeff)/sum(abs2,z.-mean(z))
            println("HISTORY ",name," start=",start," peak_fit_sigma=",coeff[2]," r2=",r2)
        end
    end
end
main()
