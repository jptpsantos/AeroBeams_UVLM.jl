The standalone LaTeX proposal is [UVLM_AEROBEAMS_MIGRATION_PROPOSAL.tex](UVLM_AEROBEAMS_MIGRATION_PROPOSAL.tex). It expands the [Markdown roadmap](UVLM_AEROBEAMS_MIGRATION.md) with explanations, terminology comparisons, mathematical distinctions, proposed API examples, and an explicit approval checkpoint.

The migration is proposed work and has not been started. The document distinguishes P0–P5 (terminology and organization) from P6–P9 (AeroBeams coupling), with performance and advanced capabilities covered separately.

To compile, open a terminal in this directory and run:

```powershell
pdflatex -interaction=nonstopmode -halt-on-error UVLM_AEROBEAMS_MIGRATION_PROPOSAL.tex
pdflatex -interaction=nonstopmode -halt-on-error UVLM_AEROBEAMS_MIGRATION_PROPOSAL.tex
```

Alternatively, use `latexmk -pdf UVLM_AEROBEAMS_MIGRATION_PROPOSAL.tex`, or upload the single `.tex` file to Overleaf and select **pdfLaTeX**. The second compiler pass resolves the table of contents and cross-references. No external figures, bibliography files, shell escape, or Julia execution are required.

No LaTeX compiler was found in the authoring environment. The source received static structure and reference checks; PDF compilation and visual layout verification remain unperformed.
