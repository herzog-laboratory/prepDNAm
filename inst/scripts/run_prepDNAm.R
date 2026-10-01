#!/usr/bin/env Rscript
#
# Run prepDNAm non-interactively from a YAML configuration file, e.g. under
# Slurm / Nextflow / Snakemake:
#
#   Rscript run_prepDNAm.R --config study1.yml [--force] [--quiet]
#
# The YAML keys are the arguments of prepDNAm::qc_config(), for example:
#
#   input: /data/idats/study1
#   run_dir: /data/runs/study1
#   array: EPICv2
#   pheno: /data/pheno/study1.csv
#   batch: chip
#   cores: 8
#   export_csv: true      # also write CSVs in the biolearn standard format
#
# Re-running on an existing run directory resumes after the last completed
# stage. Exit status is non-zero if a stage failed (see FAILED_<stage>.txt).

args <- commandArgs(trailingOnly = TRUE)
usage <- function() {
  cat("Usage: Rscript run_prepDNAm.R --config <file.yml> [--force] [--quiet]\n")
  quit(status = 2)
}
if (length(args) == 0 || any(args %in% c("-h", "--help"))) usage()
config_file <- args[match("--config", args) + 1]
if (is.na(config_file) || !file.exists(config_file)) usage()
force <- "--force" %in% args
if ("--quiet" %in% args) options(prepDNAm.quiet = TRUE)

suppressPackageStartupMessages(library(prepDNAm))
params <- yaml::read_yaml(config_file)
cfg <- do.call(qc_config, params)
res <- tryCatch(run_prepDNAm(cfg, force = force), error = function(e) {
  message(conditionMessage(e))
  quit(status = 1)
})
print(res)
