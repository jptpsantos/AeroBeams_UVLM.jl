# Chang wing-propeller lumped-inertia remesh audit

Generated: 2026-09-02 09:30:26

## Implemented interpretation

The spreadsheet stiffness is treated as a distributed spanwise constitutive field. The 16 mass/CG/inertia records, plus the zero root record, define nodal reference values. Each complete 6-by-6 nodal spatial-inertia block is divided by its source tributary length, the resulting density is linearly interpolated, and that density is integrated exactly over every target element. Each target integral is lumped at its outer node. No free structural node is massless.

The integral of each source-node hat function equals its tributary length. Therefore the transfer uses nonnegative weights and conserves the sum of the complete spatial-inertia matrix exactly. It conserves mass, first moments, and beam-axis inertia simultaneously, without separately rescaling signed components.

Workbook connectivity is important: inertia elements 52-67 connect each non-root flexible wing station to a short offset node at the same span coordinate. They are rigid mass elements, not the flexible spanwise beam elements 20-35. Thus source_16 is the literal concentrated-inertia reference. The uniform cases intentionally construct an equivalent continuous distribution from those nodal reference records.

## Mass-matrix verification

| Case | Elements | Mass (kg) | Spatial-inertia error | min eig(M) | min nodal-block eig |
|---|---:|---:|---:|---:|---:|
| source_16 | 16 | 177.32851668 | 0.000e+00 | +1.074102e-01 | +1.074102e-01 |
| uniform_16 | 16 | 177.32851668 | 1.603e-16 | +1.022110e-01 | +1.022110e-01 |
| uniform_20 | 20 | 177.32851668 | 1.603e-16 | +7.840847e-02 | +7.840847e-02 |
| uniform_24 | 24 | 177.32851668 | 1.603e-16 | +6.345043e-02 | +6.345043e-02 |
| uniform_30 | 30 | 177.32851668 | 1.603e-16 | +4.923315e-02 | +4.923315e-02 |
| uniform_32 | 32 | 177.32851668 | 1.202e-16 | +4.579587e-02 | +4.579587e-02 |
| uniform_36 | 36 | 177.32851668 | 8.014e-17 | +4.017197e-02 | +4.017197e-02 |
| uniform_40 | 40 | 177.32851668 | 4.007e-17 | +3.576775e-02 | +3.576775e-02 |
| uniform_48 | 48 | 177.32851668 | 3.206e-16 | +2.932059e-02 | +2.932059e-02 |
| uniform_60 | 60 | 177.32851668 | 1.202e-16 | +2.306575e-02 | +2.306575e-02 |
| uniform_80 | 80 | 177.32851668 | 4.808e-16 | +1.700473e-02 | +1.700473e-02 |
| uniform_120 | 120 | 177.32851668 | 1.603e-16 | +1.113905e-02 | +1.113905e-02 |
| uniform_160 | 160 | 177.32851668 | 1.202e-16 | +8.279953e-03 | +8.279953e-03 |
| uniform_240 | 240 | 177.32851668 | 4.408e-16 | +5.166785e-03 | +5.166785e-03 |
| uniform_320 | 320 | 177.32851668 | 2.003e-16 | +2.906316e-03 | +2.906316e-03 |

## Thirty-element remap interpretation sensitivity

| Method | Mass (kg) | Spatial-inertia error | min eig(M) | OOP1 | IP1 | OOP2 | T1 | OOP3 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| legacy_linear | 175.41043878 | 1.082e-02 | +4.296399e-02 | 6.63808 | 20.81901 | 36.45316 | 42.89074 | 93.40783 |
| conservative_spatial | 177.32851668 | 1.603e-16 | +4.923315e-02 | 6.63049 | 20.79605 | 36.47010 | 42.84774 | 93.52212 |
| element_overlap | 177.32851668 | 8.014e-17 | +4.296399e-02 | 6.97171 | 21.85155 | 37.48588 | 44.24148 | 94.78981 |

legacy_linear is the former point-sampled nodal-density interpolation. conservative_spatial is the corrected exact-integration version. element_overlap assumes the records belong to the preceding spanwise flexible elements; it is included only to quantify that alternative interpretation, which conflicts with the workbook connectivity.

## Isolated-wing natural frequencies

| Case | OOP1 | IP1 | OOP2 | T1 | OOP3 |
|---|---:|---:|---:|---:|---:|
| source_16 | 6.62627 | 20.78295 | 36.38983 | 42.87465 | 92.42593 |
| uniform_16 | 6.63403 | 20.80353 | 36.40711 | 42.87282 | 92.79187 |
| uniform_20 | 6.63615 | 20.81031 | 36.46820 | 42.90134 | 93.25806 |
| uniform_24 | 6.63614 | 20.81087 | 36.48578 | 42.89387 | 93.40853 |
| uniform_30 | 6.63049 | 20.79605 | 36.47010 | 42.84774 | 93.52212 |
| uniform_32 | 6.63652 | 20.81284 | 36.51431 | 42.90502 | 93.61634 |
| uniform_36 | 6.63668 | 20.81322 | 36.50816 | 42.88693 | 93.63485 |
| uniform_40 | 6.63621 | 20.81251 | 36.52306 | 42.90235 | 93.69693 |
| uniform_48 | 6.63571 | 20.81127 | 36.52234 | 42.89208 | 93.72314 |
| uniform_60 | 6.63299 | 20.80378 | 36.50564 | 42.86372 | 93.73117 |
| uniform_80 | 6.63582 | 20.81170 | 36.52569 | 42.88710 | 93.78645 |
| uniform_120 | 6.63362 | 20.80571 | 36.51434 | 42.86766 | 93.78048 |
| uniform_160 | 6.63476 | 20.80894 | 36.52376 | 42.87981 | 93.80210 |
| uniform_240 | 6.63376 | 20.80619 | 36.51652 | 42.86865 | 93.79291 |
| uniform_320 | 6.63440 | 20.80793 | 36.52059 | 42.87411 | 93.80168 |
| Published reference | 6.63000 | 20.93000 | 36.66000 | 43.50000 | 93.31000 |

The exact 16-record/source-station model differs from the published values by OOP1 -0.056%, IP1 -0.703%, OOP2 -0.737%, T1 -1.438%, OOP3 -0.947%.

The valid 20-element uniform remesh is not yet frequency-converged: its differences from the published values are OOP1 +0.093%, IP1 -0.572%, OOP2 -0.523%, T1 -1.376%, OOP3 -0.056%.

At the finest uniform mesh (320 elements), the differences from the published isolated-wing values are OOP1 +0.066%, IP1 -0.583%, OOP2 -0.380%, T1 -1.439%, OOP3 +0.527%.

## Coupled wing-pylon-propeller natural frequencies

The table gives the first 10 undamped, non-gyroscopic structural frequencies. The pitch/yaw attachment is kept at the exact span fraction 0.83 by work-conjugate interpolation between adjacent beam nodes.
This intentionally removes attachment-node snapping from the modal remesh comparison. The unmodified production driver still selects the nearest structural node; in the 20-element end-to-end test it attached at y = 6.375 m instead of the exact 6.225 m position, so that separate source of mesh dependence remains outside the inertia-remap correction.

| Case | f1 | f2 | f3 | f4 | f5 | f6 | f7 | f8 | f9 | f10 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| source_16 | 5.90034 | 7.95783 | 8.52744 | 19.38857 | 35.10018 | 39.27353 | 91.01576 | 108.18211 | 114.09417 | 155.65989 |
| uniform_16 | 5.90702 | 7.95799 | 8.52998 | 19.40898 | 35.17572 | 39.22627 | 91.47607 | 108.08763 | 114.43300 | 156.63029 |
| uniform_20 | 5.90794 | 7.95797 | 8.53077 | 19.41375 | 35.21309 | 39.26788 | 91.92720 | 108.25042 | 114.60443 | 156.97904 |
| uniform_24 | 5.90739 | 7.95785 | 8.52991 | 19.41364 | 35.20129 | 39.28422 | 91.98842 | 108.33478 | 114.54513 | 156.90473 |
| uniform_30 | 5.90287 | 7.95783 | 8.52882 | 19.39934 | 35.15803 | 39.27965 | 92.09443 | 108.36857 | 114.52226 | 156.74651 |
| uniform_32 | 5.90745 | 7.95790 | 8.53082 | 19.41465 | 35.22001 | 39.30616 | 92.24298 | 108.41373 | 114.67195 | 157.01078 |
| uniform_36 | 5.90733 | 7.95784 | 8.53032 | 19.41485 | 35.19891 | 39.30550 | 92.21677 | 108.43429 | 114.61001 | 156.91565 |
| uniform_40 | 5.90698 | 7.95785 | 8.53044 | 19.41404 | 35.21645 | 39.31173 | 92.29830 | 108.44674 | 114.65119 | 156.99044 |
| uniform_48 | 5.90648 | 7.95784 | 8.53021 | 19.41271 | 35.20804 | 39.31349 | 92.31004 | 108.46650 | 114.65242 | 156.93601 |
| uniform_60 | 5.90434 | 7.95783 | 8.52962 | 19.40563 | 35.17760 | 39.30770 | 92.30770 | 108.46869 | 114.60917 | 156.81367 |
| uniform_80 | 5.90641 | 7.95784 | 8.53034 | 19.41282 | 35.20112 | 39.31951 | 92.37313 | 108.49307 | 114.64941 | 156.91435 |
| uniform_120 | 5.90469 | 7.95782 | 8.52974 | 19.40720 | 35.18167 | 39.31392 | 92.35424 | 108.49156 | 114.61461 | 156.82598 |
| uniform_160 | 5.90551 | 7.95781 | 8.52989 | 19.41016 | 35.19319 | 39.31810 | 92.37469 | 108.49952 | 114.62183 | 156.87576 |
| uniform_240 | 5.90474 | 7.95779 | 8.52963 | 19.40760 | 35.18154 | 39.31402 | 92.36005 | 108.49504 | 114.58878 | 156.82683 |
| uniform_320 | 5.90521 | 7.95779 | 8.52978 | 19.40920 | 35.18707 | 39.31663 | 92.37101 | 108.49956 | 114.60404 | 156.85025 |

Relative to the exact 16-record/source-station model, the finest uniform-mesh changes in the first ten coupled frequencies are: f1 +0.083%, f2 -0.000%, f3 +0.027%, f4 +0.106%, f5 +0.248%, f6 +0.110%, f7 +1.489%, f8 +0.293%, f9 +0.447%, f10 +0.765%.

## Structural mesh convergence relative to the finest uniform mesh

| Elements | Maximum wing-mode error | Maximum coupled-mode error |
|---:|---:|---:|
| 16 | 1.0765% | 0.9689% |
| 20 | 0.5795% | 0.4805% |
| 24 | 0.4191% | 0.4142% |
| 30 | 0.2980% | 0.2994% |
| 32 | 0.1976% | 0.1386% |
| 36 | 0.1779% | 0.1670% |
| 40 | 0.1117% | 0.0894% |
| 48 | 0.0837% | 0.0660% |
| 60 | 0.0752% | 0.0685% |
| 80 | 0.0303% | 0.0409% |
| 120 | 0.0226% | 0.0182% |
| 160 | 0.0133% | 0.0174% |
| 240 | 0.0127% | 0.0157% |
| 320 | 0.0000% | 0.0000% |

## Why the present 20-element approach fails

| Remap | min eig(M) | min eig(recovered J_CG) | max recovered CG offset (m) | Modal solve |
|---|---:|---:|---:|---|
| legacy_linear | +6.444598e-02 | +0.000000e+00 | 0.269021 | success |
| signed_scaled | -1.874893e+01 | -5.546860e+01 | 4.146195 | PosDefException |

The signed-scaled method globally rescales each first moment and product of inertia by its own ratio of source total to interpolated total. Several of these totals are small because positive and negative entries cancel. At 20 elements, the independent scale factors amplify interpolation error, produce unrealistic CG offsets, and make recovered center-of-mass inertia tensors and the global mass matrix indefinite. The failure therefore occurs before Generalized-alpha time integration; it is not caused by the GA parameters.

The old unscaled componentwise interpolation can remain positive definite, but it does not conserve the complete spatial inertia. The conservative spatial-block remap satisfies both requirements.

## End-to-end 20-element aeroelastic validation

A separate in-memory launcher injected the conservative arrays into the actual production UVLM/Generalized-alpha driver. The short trim-plus-impulse test is a solver smoke test, not a damping-convergence result.

- integrated steps: 17
- finite history: true
- all coupling steps converged: true
- maximum coupling iterations: 2
- maximum state residual: 1.7031887696489772e-8
- maximum load residual: 1.453498553191216e-5
- maximum coupled-equilibrium residual: 2.2689832797192592e-5
- validation passed: true
