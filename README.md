
# prepDNAm

<!-- badges: start -->
<!-- badges: end -->

prepDNAm takes Illumina DNA methylation array data from raw IDAT files,
signal-intensity tables (e.g. from GEO) or beta matrices to an analysis-ready
beta matrix: sample QC, normalisation (Noob + BMIQ), probe filtering and
imputation, with an interactive HTML QC report.

- **Arrays:** 450k, EPIC v1, EPIC v2, Mouse Methylation.
- **Resumable:** every stage writes its results into a run directory and is
  skipped if they already exist, so a crashed or interrupted run continues
  where it stopped.
- **Batch-wise:** reading and normalisation run per batch (plate, chip, or any
  sample sheet column); filtering and imputation on the pooled data, so the
  result does not depend on the batching.
- **Reproducible:** the run directory records the configuration, sample sheet,
  package versions, seed, input checksums and a log.
- **Scriptable:** run from R or from the command line (HPC, Slurm, Snakemake).

## Installation

prepDNAm needs R >= 4.1 and several Bioconductor packages.

``` r
install.packages(c("remotes", "BiocManager"))

# core Bioconductor dependencies
BiocManager::install(c("minfi", "ChAMP", "impute"))

# the manifest and annotation package for your array, e.g. EPIC v2:
BiocManager::install(c("IlluminaHumanMethylationEPICv2manifest",
                       "IlluminaHumanMethylationEPICv2anno.20a1.hg38"))

remotes::install_github("herzog-laboratory/prepDNAm")
```

See `vignette("prepDNAm")` for the packages needed per array and the optional
ones (genotype checks, cell-type and age estimates, UMAP, the report).

## Quick start

``` r
library(prepDNAm)

cfg <- qc_config(
  input   = "data/idats",          # directory of IDATs (searched recursively)
  run_dir = "runs/study1",         # everything is written here
  array   = "EPICv2",
  pheno   = "data/samples.csv",    # optional sample sheet, see ?sample_sheet
  cores   = 4
)
run_prepDNAm(cfg)                  # all stages; re-run to resume
qc_report("runs/study1", open = TRUE)

beta <- readRDS("runs/study1/beta_merged.rds")
```

From the command line, with the same options in a YAML file:

``` sh
Rscript $(Rscript -e 'cat(system.file("scripts/run_prepDNAm.R", package = "prepDNAm"))') \
  --config study1.yml
```

## Documentation

- `vignette("prepDNAm")`: installation, running, outputs, the report,
  exporting, running on a cluster.
- `vignette("sample-sheet")`: how to prepare the sample sheet (sample IDs,
  sex, age, donor, batch, exclusions).
- `?qc_config`: all options and their defaults.
