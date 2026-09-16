# Implemented UVLM interface migration

Read [UVLM_ANALYSIS_GUIDE.tex](UVLM_ANALYSIS_GUIDE.tex) for the self-contained code and analysis guide. It includes both convergence stages, the single Chang case, prescribed-motion aerodynamics, the file responsibility map and embedded source listings. Compile twice with pdfLaTeX, or upload the single file to Overleaf.

The implementation introduces AeroBeams-style constructor/problem terminology over the retained numerical backend and centralizes the studies. It does not implement the later native AeroBeams structural bridge. Remaining roadmap work is stated explicitly in the guide.

[baseline.toml](baseline.toml) records the original numerical source identity. [VALIDATION.md](VALIDATION.md) records verification and limits. Existing aerodynamic selection seals are preserved; migrated runs must establish evidence under the new source identity.
