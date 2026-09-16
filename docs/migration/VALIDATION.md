# Migration verification — 15 September 2026

The final declared studies environment passed **508 assertions** with Julia 1.11.5. The runner sets `LOAD_PATH=["@","@stdlib"]`, so package resolution does not fall back to the global environment. Installed package source caches were shared through Julia's normal depot mechanism.

| Test group | Passed |
|---|---:|
| Canonical settings, entry inclusion and aerodynamic dry run | 12 |
| Existing Chang configuration and propagation | 73 |
| New model/problem API and accepted wake history | 49 |
| Existing aerodynamic library tests | 130 |
| Existing aerodynamic convergence integrity | 41 |
| Existing fixed-core study | 66 |
| Independent two-stage convergence | 37 |
| Independent aeroelastic configuration | 24 |
| Moving-block case diagnostics | 8 |
| Reference moving-block FFT method | 55 |
| Executable examples, Chang parity and premature-stop rejection | 13 |
| **Total** | **508** |

The pre-edit suite passed 434 assertions. The added coverage is 74 assertions. The no-speed-override fixture was made explicit instead of depending on the user's editable aeroelastic preset; the current 80 m/s setting is retained.

Numerical evidence includes circulation, temporal rates and force comparisons against the existing one-call propagation; wake geometry through mature-wake convection and capacity truncation; A–B–A trial determinism; rollback; duplicate-commit rejection; saved-result isolation; and steady-backend parity. A 38-step Chang response at 85 m/s produces exactly matching displacement, velocity and coupling-iteration histories through the old and new interfaces. A deliberately early-stopped study case is rejected before damping acceptance and records incomplete status.

The existing two-stage tests use synthetic, hash-verified histories to test selection integrity and handoff. They do not assert that a production aerodynamic sweep has converged. Moving-block regressions compare against an independent transcription of the selected reference method.

Environment checks:

- The aerodynamic package loads in isolation with FFTW and Plots absent from its active package load path.
- `Pkg.develop` and `Pkg.resolve` succeeded offline using the installed registry and package cache. The machine-generated studies manifest records a relative `..` path to WingPropellerUVLM.
- The lockfile was checked against both packages' declared compatibility bounds, including Interpolations 0.15.1.
- Julia's compressed-registry subprocess initially hit a Windows `EBADF` error. Unpacking the installed registry into a temporary depot allowed normal package resolution to complete. Temporary-depot pidfile cleanup emitted Windows `EINVAL` warnings; resolution and tests exited successfully. The user's global project and registry were not modified.

Document checks:

- The self-contained guide embeds ten complete entry/interface source listings; each was compared with its actual file.
- LaTeX brace balance, environment counts and local document links were checked.
- No pdfLaTeX, latexmk or tectonic executable was available. A PDF was not compiled or visually inspected.
- `git diff --check` passed.

Reproduce the full suite from the repository root after study setup:

```sh
julia --startup-file=no --project=lib/WingPropellerUVLM/studies lib/WingPropellerUVLM/studies/test/runtests.jl
```

The authoring run used `--compiled-modules=existing` and a temporary first depot to avoid writes to the restricted global depot. [test_output.log](test_output.log) preserves the successful run; [verification.toml](verification.toml) records the final source and estimator hashes.

No full six-second damping convergence sweep, production aerodynamic selection, native AeroBeams structural coupling, or PDF rendering is claimed. Original selections/seals remain unchanged. Remaining roadmap work is listed explicitly in [UVLM_ANALYSIS_GUIDE.tex](UVLM_ANALYSIS_GUIDE.tex).
