# Chang reference data and model construction. Helpers below are pure numerical
# operations; build_chang_model_parameters(config) owns all derived case data.

chang_skew3(v) = [
    0.0 -v[3] v[2];
    v[3] 0.0 -v[1];
    -v[2] v[1] 0.0
]

"""Build a complete nodal spatial inertia about the beam reference axis."""
function chang_spatial_inertia_block(mass,inertia_at_cg,offset)
    skew_offset = chang_skew3(offset)
    identity3 = Matrix{Float64}(I,3,3)
    rotationalInertia = inertia_at_cg .- mass .* skew_offset * skew_offset
    return [mass.*identity3 -mass.*skew_offset;
        mass.*skew_offset rotationalInertia]
end

"""
    chang_control_volume_spatial_remap(source_nodes, source_blocks, target_nodes)

Convert the source nodal spatial inertias to a piecewise-linear matrix density
and integrate it exactly over target nodal control volumes.  Complete 6-by-6
blocks are transferred with nonnegative weights, so mass, first moments, all
products of inertia, symmetry, and positive semidefiniteness remain coupled.
The root control-volume contribution is moved to the first free node because
the root DOFs are constrained and the reference root inertia is zero.
"""
function chang_control_volume_spatial_remap(
    source_nodes::AbstractVector{<:Real},
    source_blocks::AbstractVector{<:AbstractMatrix},
    target_nodes::AbstractVector{<:Real}
)
    length(source_nodes) == length(source_blocks) || throw(DimensionMismatch(
        "source_nodes and source_blocks must have the same length"
    ))
    issorted(source_nodes) || throw(ArgumentError("source_nodes must be sorted"))
    issorted(target_nodes) || throw(ArgumentError("target_nodes must be sorted"))
    isapprox(target_nodes[1],source_nodes[1]; atol=1e-12,rtol=0.0) ||
        throw(ArgumentError("source and target meshes must have the same root"))
    isapprox(target_nodes[end],source_nodes[end]; atol=1e-10,rtol=0.0) ||
        throw(ArgumentError("source and target meshes must have the same tip"))

    sourceLengths = zeros(length(source_nodes))
    sourceElementLengths = diff(source_nodes)
    sourceLengths[1] = sourceElementLengths[1] / 2
    sourceLengths[end] = sourceElementLengths[end] / 2
    for node = 2:(length(source_nodes)-1)
        sourceLengths[node] =
            (sourceElementLengths[node-1] + sourceElementLengths[node]) / 2
    end
    sourceDensity = [
        Matrix{Float64}(source_blocks[node]) ./ sourceLengths[node]
        for node in eachindex(source_nodes)
    ]

    controlEdges = [
        target_nodes[1];
        (target_nodes[1:end-1] .+ target_nodes[2:end]) ./ 2;
        target_nodes[end]
    ]
    targetBlocks = [zeros(6,6) for _ in eachindex(target_nodes)]

    for targetNode in eachindex(target_nodes)
        targetLeft = controlEdges[targetNode]
        targetRight = controlEdges[targetNode+1]
        targetBlock = targetBlocks[targetNode]

        for sourceElement = 1:(length(source_nodes)-1)
            sourceLeft = source_nodes[sourceElement]
            sourceRight = source_nodes[sourceElement+1]
            overlapLeft = max(targetLeft,sourceLeft)
            overlapRight = min(targetRight,sourceRight)
            overlapRight <= overlapLeft && continue

            interval = sourceRight - sourceLeft
            tLeft = (overlapLeft - sourceLeft) / interval
            tRight = (overlapRight - sourceLeft) / interval
            leftWeight = interval * (
                (tRight - tRight^2 / 2) - (tLeft - tLeft^2 / 2)
            )
            rightWeight = interval * (tRight^2 - tLeft^2) / 2
            targetBlock .+= leftWeight .* sourceDensity[sourceElement]
            targetBlock .+= rightWeight .* sourceDensity[sourceElement+1]
        end
        targetBlock .= 0.5 .* (targetBlock .+ targetBlock')
    end

    # The clamped root block would be discarded during boundary-condition
    # elimination. Move it to the first free station to retain the complete
    # reconstructed wing inertia in the dynamic model.
    targetBlocks[2] .+= targetBlocks[1]
    fill!(targetBlocks[1],0.0)
    return targetBlocks
end

function chang_properties_from_spatial_inertia(block)
    mass = tr(block[1:3,1:3]) / 3
    mass <= 100eps(Float64) && return (
        mass=0.0,
        cg=zeros(3),
        inertia_at_cg=zeros(3,3)
    )
    skew_offset = -block[1:3,4:6] ./ mass
    skew_offset = 0.5 .* (skew_offset .- skew_offset')
    offset = [skew_offset[3,2],skew_offset[1,3],skew_offset[2,1]]
    inertia_at_cg = block[4:6,4:6] .+ mass .* skew_offset * skew_offset
    inertia_at_cg = 0.5 .* (inertia_at_cg .+ inertia_at_cg')
    return (; mass,cg=offset,inertia_at_cg)
end

"""Build geometry, reference distributions, rotor properties, and time grid from one case."""
function build_chang_model_parameters(config)
    # Chang-specific geometry, structural distributions, propeller properties,
    # and simulation environment. Reusable algorithms belong in WingPropellerUVLM.
    # Run settings come from chang_case.jl; the tables below are Chang reference data.

    # ### Wing model
    Cr = config.wing.root_chord_m
    Ct = config.wing.tip_chord_m
    b = config.wing.span_m
    span_length = b
    taperRatio = Ct / Cr
    S_wing = b * (Cr + Ct) / 2.0
    xle_tip = 0.50 * (Cr - Ct)
    xle = [0.0,xle_tip]
    yle = [0.0,b]
    zle = [0.0,0.0]
    chord_geo = [Cr,Ct]
    theta_geo = [0,0]
    phi_geo = [0,0]
    spacing_s_wing = Uniform()
    spacing_c_wing = Uniform()
    mirror_wing = false
    symmetric_wing = config.wing.symmetric

    ns_wing = config.wing.spanwise_panels
    nc_wing = config.wing.chordwise_panels

    Ne = ns_wing
    le = span_length / Ne
    nnodes = Ne + 1
    ndof = 6
    NDOF = nnodes * ndof

    X_nodes_new = collect(range(0,span_length,length=nnodes))
    X_cen_new = [(i - 0.5) * le for i = 1:Ne]
    etaCenNew = X_cen_new ./ span_length

    chord = collect(range(Cr,Ct,length=nnodes))
    xle_distribution = collect(range(xle[1],xle[2],length=nnodes))

    span_nodes = X_nodes_new
    L_node_new = fill(le,nnodes)
    L_node_new[1] = le / 2.0
    L_node_new[end] = le / 2.0

    # --- 3.2 Nodal Mass Data Processing ---
    changFtToM = 0.3048
    changSlugToKg = 14.5939029372064
    changSlugFt2ToKgm2 = changSlugToKg * changFtToM^2
    chang_lbf_to_n = 4.4482216152605
    changLbfFt2ToNm2 = chang_lbf_to_n * changFtToM^2

    X_nodes_ft = [
        0.0,
        0.54675,
        1.0936,
        2.18722,
        3.28084,
        4.51115,
        5.74147,
        6.97178,
        8.2021,
        10.2526,
        12.30315,
        14.353675,
        16.4042,
        18.4547,
        20.50525,
        22.55577,
        24.6063
    ]
    X_nodes_orig_unscaled = X_nodes_ft .* changFtToM
    # The workbook tip coordinate differs from 7.5 m only by printed precision.
    # Scaling all source stations by the same factor makes the remap cover exactly
    # the active [0, span_length] interval and also keeps the data usable if the
    # reference span is changed parametrically.
    X_nodes_orig = X_nodes_orig_unscaled .* (span_length / X_nodes_orig_unscaled[end])
    L_orig_elements = diff(X_nodes_orig)

    L_node_orig = zeros(17)
    L_node_orig[1] = L_orig_elements[1] / 2.0
    for i = 2:16
        L_node_orig[i] = (L_orig_elements[i-1] + L_orig_elements[i]) / 2.0
    end
    L_node_orig[17] = L_orig_elements[16] / 2.0

    Mass_slug = [
        0.4276211,
        0.4139885,
        0.736678,
        0.4139886,
        0.9436633,
        0.7336013,
        0.3851399,
        0.9932094,
        0.9932177,
        0.9932237,
        0.9931985,
        0.971823,
        0.9687671,
        0.9109462,
        0.9367904,
        0.3350064
    ]
    Ixx_slug = [
        0.6055,
        0.6281,
        1.92,
        0.6281,
        2.3536,
        1.909,
        0.5645,
        2.4813,
        2.4813,
        2.4814,
        2.4813,
        2.3998,
        2.3877,
        2.2748,
        2.3653,
        0.459
    ]
    # Keep the inertia components in the literal Model Properties.xlsx column
    # order: X is spanwise, Y is chordwise, and Z is vertical.  The later
    # structural-to-aerodynamic adapter performs the required basis conversion;
    # applying another Y/Z permutation here is incorrect.
    Iyy_slug = [
        0.0969,
        0.1028,
        0.1269,
        0.1028,
        0.2878,
        0.1266,
        0.0967,
        0.3565,
        0.3656,
        0.3656,
        0.3656,
        0.3577,
        0.3561,
        0.3028,
        0.3572,
        0.0797
    ]
    Izz_slug = [
        0.6475,
        0.6777,
        1.9703,
        0.6777,
        2.3194,
        1.9593,
        0.6118,
        2.388,
        2.388,
        2.388,
        2.3879,
        2.3118,
        2.3009,
        2.2058,
        2.2509,
        0.4844
    ]
    Ixy_slug = [
        0.000483499,
        0.000545311,
        -0.003364902,
        0.000546007,
        -0.003767303,
        -0.003368369,
        0.000500823,
        -0.004768277,
        -0.004781691,
        -0.004764693,
        -0.00476285,
        -0.004980494,
        -0.004983901,
        -0.005318305,
        -0.005590199,
        0.000500334
    ]
    Ixz_slug = [
        0.001076907,
        -1.83141E-07,
        0.001721633,
        -4.2213E-07,
        0.004030091,
        -0.001715619,
        0.000169958,
        0.002555006,
        0.002553573,
        0.002555579,
        0.002557012,
        0.002683975,
        0.00262101,
        0.003720733,
        0.004553367,
        -0.001762152
    ]
    Iyz_slug = [
        -0.000821912,
        6.13447E-06,
        0.026378458,
        7.88634E-06,
        0.061766804,
        -0.026536753,
        -7.79545E-06,
        0.031499987,
        0.031497299,
        0.03149642,
        0.031500523,
        0.030281192,
        0.03089678,
        0.024007273,
        0.035641812,
        0.001067331
    ]

    cgXFt = [
        0.027173898,
        0.037438118,
        0.021032451,
        0.037438576,
        0.02464125,
        0.021088361,
        0.040243862,
        0.031203696,
        0.03121296,
        0.03120211,
        0.031192549,
        0.031900145,
        0.031983194,
        0.034021895,
        0.03582014,
        0.039275863
    ]
    cgYFt = [
        -0.66071232,
        -0.668605558,
        -0.919885764,
        -0.668609685,
        -0.864794967,
        -0.921385197,
        -0.671753688,
        -0.857220095,
        -0.857206007,
        -0.857200098,
        -0.857224387,
        -0.863614669,
        -0.864554528,
        -0.874965167,
        -0.869545335,
        -0.664042589
    ]
    cgZFt = [
        -0.245439677,
        -0.17398455,
        0.016601162,
        -0.014913878,
        0.121183138,
        -0.199929292,
        0.137702621,
        -0.081678808,
        -0.059307771,
        -0.03693769,
        -0.014573273,
        0.01200473,
        0.03224826,
        0.09014999,
        0.048423281,
        0.248940617
    ]

    M_orig_nodes = [0.0; Mass_slug .* changSlugToKg]
    Ixx_orig_nodes = [0.0; Ixx_slug .* changSlugFt2ToKgm2]
    Iyy_orig_nodes = [0.0; Iyy_slug .* changSlugFt2ToKgm2]
    Izz_orig_nodes = [0.0; Izz_slug .* changSlugFt2ToKgm2]
    Ixy_orig_nodes = [0.0; Ixy_slug .* changSlugFt2ToKgm2]
    Ixz_orig_nodes = [0.0; Ixz_slug .* changSlugFt2ToKgm2]
    Iyz_orig_nodes = [0.0; Iyz_slug .* changSlugFt2ToKgm2]

    cg_x_orig_nodes = [cgXFt[1] * changFtToM; cgXFt .* changFtToM]
    cg_y_orig_nodes = [cgYFt[1] * changFtToM; cgYFt .* changFtToM]
    cg_z_orig_nodes = [cgZFt[1] * changFtToM; cgZFt .* changFtToM]

    source_spatial_inertia_blocks = Matrix{Float64}[]
    for node in eachindex(X_nodes_orig)
        inertia_at_cg = [
            Ixx_orig_nodes[node] Ixy_orig_nodes[node] Ixz_orig_nodes[node];
            Ixy_orig_nodes[node] Iyy_orig_nodes[node] Iyz_orig_nodes[node];
            Ixz_orig_nodes[node] Iyz_orig_nodes[node] Izz_orig_nodes[node]
        ]
        offset = [
            cg_x_orig_nodes[node],
            cg_y_orig_nodes[node],
            cg_z_orig_nodes[node]
        ]
        push!(
            source_spatial_inertia_blocks,
            chang_spatial_inertia_block(M_orig_nodes[node],inertia_at_cg,offset)
        )
    end

    spatial_inertia_node_blocks = chang_control_volume_spatial_remap(
        X_nodes_orig,
        source_spatial_inertia_blocks,
        X_nodes_new
    )

    remappedInertiaProperties =
        chang_properties_from_spatial_inertia.(spatial_inertia_node_blocks)
    m_node_vec = [property.mass for property in remappedInertiaProperties]
    cg_x_node_vec = [property.cg[1] for property in remappedInertiaProperties]
    cg_y_node_vec = [property.cg[2] for property in remappedInertiaProperties]
    cg_z_node_vec = [property.cg[3] for property in remappedInertiaProperties]
    Ixx_node_vec = [property.inertia_at_cg[1,1] for property in remappedInertiaProperties]
    Iyy_node_vec = [property.inertia_at_cg[2,2] for property in remappedInertiaProperties]
    Izz_node_vec = [property.inertia_at_cg[3,3] for property in remappedInertiaProperties]
    Ixy_node_vec = [property.inertia_at_cg[1,2] for property in remappedInertiaProperties]
    Ixz_node_vec = [property.inertia_at_cg[1,3] for property in remappedInertiaProperties]
    Iyz_node_vec = [property.inertia_at_cg[2,3] for property in remappedInertiaProperties]

    # Retain conventional first-moment and beam-axis inertia arrays for diagnostics
    # and for the explicit legacy-inertia comparison.
    mass_moment_x_node = m_node_vec .* cg_x_node_vec
    mass_moment_y_node = m_node_vec .* cg_y_node_vec
    mass_moment_z_node = m_node_vec .* cg_z_node_vec
    Ixx_beam_orig = [block[4,4] for block in source_spatial_inertia_blocks]
    Iyy_beam_orig = [block[5,5] for block in source_spatial_inertia_blocks]
    Izz_beam_orig = [block[6,6] for block in source_spatial_inertia_blocks]
    Ixy_beam_orig = [block[4,5] for block in source_spatial_inertia_blocks]
    Ixz_beam_orig = [block[4,6] for block in source_spatial_inertia_blocks]
    Iyz_beam_orig = [block[5,6] for block in source_spatial_inertia_blocks]
    Ixx_beam_node = [block[4,4] for block in spatial_inertia_node_blocks]
    Iyy_beam_node = [block[5,5] for block in spatial_inertia_node_blocks]
    Izz_beam_node = [block[6,6] for block in spatial_inertia_node_blocks]
    Ixy_beam_node = [block[4,5] for block in spatial_inertia_node_blocks]
    Ixz_beam_node = [block[4,6] for block in spatial_inertia_node_blocks]
    Iyz_beam_node = [block[5,6] for block in spatial_inertia_node_blocks]

    sourceSpatialInertiaTotal = reduce(+,source_spatial_inertia_blocks)
    remappedSpatialInertiaTotal = reduce(+,spatial_inertia_node_blocks)
    wingSpatialInertiaRelativeError =
        norm(
            remappedSpatialInertiaTotal - sourceSpatialInertiaTotal,
            Inf
        ) / max(norm(sourceSpatialInertiaTotal,Inf),1.0)
    wingSpatialInertiaRelativeError <= 1e-12 || error(
        "Wing spatial-inertia remap is not conservative: " *
        "relative error=$wingSpatialInertiaRelativeError"
    )
    all(m_node_vec[2:end] .> 0.0) || error(
        "Every free structural wing node must have a positive remapped mass"
    )
    wingSpatialBlockMinimumEigenvalue = minimum(
        minimum(eigvals(Symmetric(block)))
        for block in spatial_inertia_node_blocks[2:end]
    )
    wingSpatialBlockMinimumEigenvalue > 0.0 || error(
        "A remapped free-node spatial inertia is not positive definite: " *
        "minimum eigenvalue=$wingSpatialBlockMinimumEigenvalue"
    )
    println(
        "Wing inertia remap: control-volume spatial blocks, Ne=$Ne, " *
        "mass=$(sum(m_node_vec)) kg, relative error=$wingSpatialInertiaRelativeError, " *
        "minimum block eigenvalue=$wingSpatialBlockMinimumEigenvalue"
    )

    # --- 3.3 Full 28-Element Stiffness Data ---
    etaStiff = [
        0.04443,
        0.04444,
        0.08888,
        0.08889,
        0.13332,
        0.13333,
        0.18332,
        0.18333,
        0.23332,
        0.23333,
        0.28332,
        0.28333,
        0.33333,
        0.33334,
        0.41666,
        0.41667,
        0.49999,
        0.5,
        0.58332,
        0.58333,
        0.66666,
        0.66667,
        0.74999,
        0.75,
        0.83333,
        0.83334,
        0.91667,
        1.0
    ]

    EIyy_lbf = [
        32836890.0,
        32836890.0,
        32836890.0,
        32836890.0,
        32836890.0,
        20120070.0,
        20120070.0,
        20120070.0,
        20120070.0,
        16181100.0,
        16181100.0,
        16181100.0,
        16181100.0,
        19251360.0,
        19251360.0,
        19251360.0,
        19251360.0,
        18286380.0,
        18286380.0,
        18286380.0,
        18286380.0,
        17183860.0,
        17183860.0,
        17183860.0,
        17183860.0,
        11198940.0,
        11198940.0,
        11198940.0
    ]
    EIzz_lbf = [
        325381340.0,
        325381340.0,
        325381340.0,
        325381340.0,
        325381340.0,
        208531310.0,
        208531310.0,
        208531310.0,
        208531310.0,
        163131180.0,
        163131180.0,
        163131180.0,
        163131180.0,
        194913580.0,
        194913580.0,
        194913580.0,
        194913580.0,
        181528760.0,
        181528760.0,
        181528760.0,
        181528760.0,
        163800810.0,
        163800810.0,
        163800810.0,
        163800810.0,
        128273870.0,
        128273870.0,
        128273870.0
    ]
    EIzy_lbf = [
        4.7103163E+06,
        4.7103163E+06,
        4.7103163E+06,
        4.7103163E+06,
        4.7103163E+06,
        2.9685812E+06,
        2.9685812E+06,
        2.9685812E+06,
        2.9685812E+06,
        2.5060481E+06,
        2.5060481E+06,
        2.5060481E+06,
        2.5060481E+06,
        3.3875434E+06,
        3.3875434E+06,
        3.3875434E+06,
        3.3875434E+06,
        3.1960560E+06,
        3.1960560E+06,
        3.1960560E+06,
        3.1960560E+06,
        3.3050893E+06,
        3.3050893E+06,
        3.3050893E+06,
        3.3050893E+06,
        2.4234529E+06,
        2.4234529E+06,
        2.4234529E+06
    ]
    GJ_lbf = [
        3.05E+07,
        3.05E+07,
        3.05E+07,
        3.05E+07,
        3.05E+07,
        1.84E+07,
        1.84E+07,
        1.84E+07,
        1.84E+07,
        1.40E+07,
        1.40E+07,
        1.40E+07,
        1.40E+07,
        1.62E+07,
        1.62E+07,
        1.62E+07,
        1.62E+07,
        1.52E+07,
        1.52E+07,
        1.52E+07,
        1.52E+07,
        1.30E+07,
        1.30E+07,
        1.30E+07,
        1.30E+07,
        9.82E+06,
        9.82E+06,
        9.82E+06
    ]
    EA_lbf = [
        168148740.0,
        168148740.0,
        168148740.0,
        168148740.0,
        168148740.0,
        111286380.0,
        111286380.0,
        111286380.0,
        111286380.0,
        85691907.0,
        85691907.0,
        85691907.0,
        85691907.0,
        102776990.0,
        102776990.0,
        102776990.0,
        102776990.0,
        96393539.0,
        96393539.0,
        96393539.0,
        96393539.0,
        88126997.0,
        88126997.0,
        88126697.0,
        88126997.0,
        100925580.0,
        100925580.0,
        100925580.0
    ]

    EIy_si = EIyy_lbf .* changLbfFt2ToNm2
    EIz_si = EIzz_lbf .* changLbfFt2ToNm2
    EIzy_si = EIzy_lbf .* changLbfFt2ToNm2
    GJ_si = GJ_lbf .* changLbfFt2ToNm2
    EA_si = EA_lbf .* chang_lbf_to_n

    EIy_vec = linear_interpolate_1d(etaStiff,EIy_si,etaCenNew)
    EIz_vec = linear_interpolate_1d(etaStiff,EIz_si,etaCenNew)
    EIzy_vec = linear_interpolate_1d(etaStiff,EIzy_si,etaCenNew)
    GJ_vec = linear_interpolate_1d(etaStiff,GJ_si,etaCenNew)
    EA_vec = linear_interpolate_1d(etaStiff,EA_si,etaCenNew)

    # The element-center vectors remain available for diagnostics and legacy
    # comparisons. The active structural assembly receives the complete
    # distribution below and integrates B'C(x)B through all stiffness breakpoints.
    wing_stiffness_distribution = (
        eta=etaStiff,
        EIy=EIy_si,
        EIz=EIz_si,
        EIzy=EIzy_si,
        GJ=GJ_si,
        EA=EA_si
    )

    # ### Propeller and pylon model
    Npropellers = length(config.propeller.attachment_eta)
    Vinf = config.simulation.freestream_speed_mps  # increase to find the flutter speed
    R_prop = config.propeller.radius_m
    Ω_wind = config.propeller.rotation_rpm
    V_trim = config.propeller.trim_speed_mps
    μ_prop = V_trim / (Ω_wind * 2 * pi / 60) / (R_prop)
    J = μ_prop * pi
    Ω = Npropellers == 0 ? 0.0 : Vinf / R_prop / μ_prop
    Nb_prop = config.propeller.blades
    c_prop = config.propeller.chord_m

    fθ_prop = config.structural.pitch_frequency_hz * 2 * pi
    fψ_prop = config.structural.yaw_frequency_hz * 2 * pi
    Kθ_prop = config.structural.pitch_stiffness_nm_per_rad
    Kψ_prop = config.structural.yaw_stiffness_nm_per_rad

    fTwist = config.structural.twist_frequency_hz * 2 * pi
    K_twist = config.structural.twist_stiffness_nm_per_rad

    ξ_prop = config.structural.propeller_damping_ratio
    # Global stiffness-proportional Rayleigh damping. For C = βK, the modal
    # damping is ξ(ω) = βω/2. Use the configured ratio at the reference frequency.
    stiffness_damping_ratio = config.structural.stiffness_damping_ratio
    stiffness_damping_reference_omega =
        isnothing(config.structural.damping_reference_frequency_hz) ?
        fθ_prop : config.structural.damping_reference_frequency_hz * 2 * pi
    η0_prop = 0.0 / R_prop

    ns_prop = config.propeller.radial_panels
    nc_prop = config.propeller.chordwise_panels

    # The Chang data already define the full blade-angle distribution. Keep this
    # offset at zero unless a deliberate collective-pitch variation is required.
    collective_pitch_offset_deg = config.propeller.collective_pitch_offset_deg
    nodal_radii,_,twists_at_nodes = get_nodal_properties_chang(ns_prop)
    blade_twists_prop = 1.0 .* (
        (twists_at_nodes .- 90 .+ collective_pitch_offset_deg) .|> deg2rad
    )

    slugToKg = changSlugToKg
    ftToM = changFtToM

    L_pylon = config.structural.pylon_length_m
    mPylon = config.structural.pylon_mass_per_length_kgpm * L_pylon
    mBlade = config.structural.blade_mass_kg
    mRotor = Nb_prop * mBlade
    mP_prop = mPylon + mRotor

    cgRotor = -L_pylon
    cgPylon = -0.5 * L_pylon

    S_root = -(mRotor * L_pylon + mPylon * L_pylon / 2.0)
    S_modal = -(mRotor * L_pylon / 2.0 + mPylon * L_pylon / 6.0)
    I_root = mRotor * L_pylon^2 + mPylon * L_pylon^2 / 3.0
    I_cross = mRotor * L_pylon^2 / 2.0 + mPylon * L_pylon^2 / 8.0

    Inθ_prop = Kθ_prop / (fθ_prop^2)
    Inψ_prop = Kψ_prop / (fψ_prop^2)
    Ix_prop = K_twist / (fTwist^2)     # positive axial inertia; Cgyro carries the sign convention

    d_EA_to_pivot = 0.0

    SθP_prop = S_modal
    SψP_prop = S_modal
    SαP_prop = S_root + mP_prop * d_EA_to_pivot
    SγP_prop = S_root + mP_prop * d_EA_to_pivot
    IθαP_prop = I_cross + d_EA_to_pivot * SθP_prop
    IψγP_prop = I_cross + d_EA_to_pivot * SψP_prop
    IαP_prop = I_root + 2.0 * d_EA_to_pivot * SθP_prop + mP_prop * d_EA_to_pivot^2
    IγP_prop = I_root + 2.0 * d_EA_to_pivot * SψP_prop + mP_prop * d_EA_to_pivot^2

    if Npropellers > 0
        println("Flexible-pylon modal structural definition:")
        println("  L_pylon = $(round(L_pylon,digits=6)) m,  mP = $(round(mP_prop,digits=6)) kg")
        println("  S_root = $(round(S_root,digits=6)), S_modal = $(round(S_modal,digits=6))")
        println("  I_root = $(round(I_root,digits=6)), I_cross = $(round(I_cross,digits=6))")
        println("  Ix = $(round(Ix_prop,digits=6))  (positive axial inertia; Cgyro sign follows spin/yaw convention)")
    end

    propeller_eta = config.propeller.attachment_eta
    propeller_span_positions = propeller_eta .* span_length
    prop_attach_nodes = [clamp(round(Int,Lp / le) + 1,1,nnodes) for Lp in propeller_span_positions]
    ndof_P = 2 * Npropellers

    for ip = 1:Npropellers
        println(
            "Propeller P$(ip) at η=$(propeller_eta[ip]), y=$(propeller_span_positions[ip]) m; " *
            "structural node=$(prop_attach_nodes[ip])"
        )
    end
    structuralPivotProp = SVector(0.0,0.0,0.0)

    # ### Simulation environment and time grid
    azimutalStepDeg = config.simulation.azimuth_step_deg
    tEnd = config.simulation.end_time_s
    tEnd > 0 || error("CHANG_END_TIME_S must be positive")
    tStep = Npropellers == 0 ? config.simulation.time_step_s :
            (azimutalStepDeg * pi / 180) / Ω
    t = range(0.0,tEnd,step=tStep)
    dt = fill(tStep,length(t) - 1)
    length(t) >= 2 || error("Simulation duration must cover at least one time step")
    println(
        Npropellers == 0 ? "Isolated-wing time step: dt = $tStep s" :
        "Timestep based on propeller rotation: dt = $tStep s"
    )

    INTERACTION_ON = config.simulation.interaction_on
    INTERACTION_MATRIX = nothing
    println("Aerodynamic Interaction: $INTERACTION_ON")

    alpha_deg = config.simulation.angle_of_attack_deg
    alpha = alpha_deg * pi / 180
    beta = config.simulation.sideslip_deg * pi / 180
    Omega_fs = [0.0,0.0,0.0]
    fs = Freestream(Vinf,alpha,beta,Omega_fs)

    Sref = S_wing
    crefWing = (2 / 3) * Cr * (1 + taperRatio + taperRatio^2) / (1 + taperRatio)
    cref = crefWing
    bref = b
    rref = [0.3 * Cr,0.0,0.0]
    ref = Reference(Sref,cref,bref,rref,Vinf,config.simulation.air_density_kgpm3)

    return (;
        b,span_length,xle,yle,zle,
        chord_geo,theta_geo,phi_geo,spacing_s_wing,spacing_c_wing,
        mirror_wing,symmetric_wing,ns_wing,nc_wing,Ne,
        le,nnodes,ndof,NDOF,chord,
        xle_distribution,span_nodes,M_orig_nodes,cg_x_orig_nodes,cg_y_orig_nodes,
        cg_z_orig_nodes,spatial_inertia_node_blocks,m_node_vec,cg_x_node_vec,cg_y_node_vec,
        cg_z_node_vec,Ixx_node_vec,Iyy_node_vec,Izz_node_vec,Ixy_node_vec,
        Ixz_node_vec,Iyz_node_vec,mass_moment_x_node,mass_moment_y_node,mass_moment_z_node,
        Ixx_beam_orig,Iyy_beam_orig,Izz_beam_orig,Ixy_beam_orig,Ixz_beam_orig,
        Iyz_beam_orig,Ixx_beam_node,Iyy_beam_node,Izz_beam_node,Ixy_beam_node,
        Ixz_beam_node,Iyz_beam_node,EIy_vec,EIz_vec,EIzy_vec,
        GJ_vec,EA_vec,wing_stiffness_distribution,Vinf,R_prop,
        Ω,Nb_prop,c_prop,Kθ_prop,Kψ_prop,
        ξ_prop,stiffness_damping_ratio,stiffness_damping_reference_omega,ns_prop,nc_prop,
        blade_twists_prop,L_pylon,mP_prop,Inθ_prop,Inψ_prop,
        Ix_prop,SθP_prop,SψP_prop,SαP_prop,SγP_prop,
        IθαP_prop,IψγP_prop,IαP_prop,IγP_prop,propeller_span_positions,
        Npropellers,prop_attach_nodes,ndof_P,
        t,dt,INTERACTION_ON,alpha,fs,
        Sref,cref,ref
    )
end
