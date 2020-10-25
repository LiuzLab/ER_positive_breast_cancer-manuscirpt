# CRISPR barcoded ER+ breast cancer metastasis analysis

Code in this repository is used to reproduce the CRISPR barcoding analysis in manuscript: **Phenotypic plasticity of ER+ breast cancer in the bone microenvironment.**

## Files description

The CRISPR barcoding analysis first uses [traceQC](https://github.com/LiuzLab/TraceQC/) package for quality control and data processing. Running the code in file ```run_traceqc.R``` outputs quality control report and aligned sequence.

Use ```glasso.R``` to perform the graphoc Lasso analysis.
