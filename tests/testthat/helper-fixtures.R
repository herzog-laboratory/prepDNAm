# Test fixtures (DESIGN.md §14.1).
#
# Everything is generated from minfiData / minfiDataEPIC into a per-session
# temp directory, lazily on first use and memoised, so test files that do not
# need the heavy fixtures (signal matrix, beta matrix) pay nothing. Nothing is
# committed. Tests call skip_if_no_fixture_data() first.

.fixtures <- new.env(parent = emptyenv())

fixture_root <- function() {
  if (is.null(.fixtures$root)) {
    .fixtures$root <- withr::local_tempdir("prepDNAm-fixtures",
                                           .local_envir = teardown_env())
  }
  .fixtures$root
}

# Memoise a fixture builder under `name`; `build(dir)` writes into `dir` and
# returns what callers get back (usually a path).
fixture_cached <- function(name, build) {
  if (is.null(.fixtures[[name]])) {
    dir <- file.path(fixture_root(), name)
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    .fixtures[[name]] <- build(dir)
  }
  .fixtures[[name]]
}

skip_if_no_fixture_data <- function(epic = FALSE) {
  skip_if_not_installed("minfiData")
  skip_if_not_installed("IlluminaHumanMethylation450kmanifest")
  if (epic) {
    skip_if_not_installed("minfiDataEPIC")
    skip_if_not_installed("IlluminaHumanMethylationEPICmanifest")
  }
}

# ---- IDATs -------------------------------------------------------------------

# minfiData ships 6 x 450k samples on two chips, one subdirectory per chip.
# `layout`: "nested" (as shipped, read-only), "flat" (all IDATs in one dir),
# "gz" (flat, gzipped).
fixture_idat_450k <- function(layout = c("nested", "flat", "gz")) {
  layout <- match.arg(layout)
  src <- system.file("extdata", package = "minfiData")
  if (layout == "nested") return(src)
  fixture_cached(paste0("idat_450k_", layout), function(dir) {
    files <- list.files(src, pattern = "\\.idat$", recursive = TRUE,
                        full.names = TRUE)
    for (f in files) {
      dest <- file.path(dir, basename(f))
      if (layout == "gz") {
        con <- gzfile(paste0(dest, ".gz"), "wb")
        writeBin(readBin(f, "raw", file.size(f)), con)
        close(con)
      } else {
        file.copy(f, dest)
      }
    }
    dir
  })
}

# minfiDataEPIC ships 3 EPIC samples, one chip each (nested layout).
fixture_idat_epic <- function() {
  system.file("extdata", package = "minfiDataEPIC")
}

# Basenames of the 450k fixture samples, in sample-sheet order.
fixture_basenames_450k <- function() {
  ss <- fixture_samplesheet_450k()
  paste(ss$Sentrix_ID, ss$Sentrix_Position, sep = "_")
}

fixture_samplesheet_450k <- function() {
  utils::read.csv(file.path(fixture_idat_450k(), "SampleSheet.csv"), skip = 7,
                  stringsAsFactors = FALSE)
}

# ---- pheno files -------------------------------------------------------------

# A pheno data frame with `basename`, `plate` (= chip), `sex`, `age`, `group`.
# `variant`: "ok" (all 6 samples), "missing" (one sample listed that has no
# IDATs), "extra" (one fixture sample left out, so it is on disk but not in
# pheno), "nobasename" (basename column renamed).
fixture_pheno_df <- function(variant = c("ok", "missing", "extra", "nobasename")) {
  variant <- match.arg(variant)
  ss <- fixture_samplesheet_450k()
  df <- data.frame(
    basename = paste(ss$Sentrix_ID, ss$Sentrix_Position, sep = "_"),
    plate = as.character(ss$Sentrix_ID),
    sex = ss$sex,
    age = ss$age,
    group = ss$Sample_Group,
    person = ss$person,
    stringsAsFactors = FALSE)
  switch(variant,
    ok = df,
    missing = rbind(df, data.frame(basename = "5723646053_R01C01", plate = "5723646053",
                                   sex = "F", age = 40, group = "GroupA", person = "id9")),
    extra = df[-nrow(df), ],
    nobasename = { names(df)[1] <- "Sample_ID"; df })
}

# Write the pheno data frame in one of the supported formats; returns the path.
fixture_pheno <- function(format = c("csv", "tsv", "txt", "xlsx", "rds", "Rdata"),
                          variant = "ok") {
  format <- match.arg(format)
  df <- fixture_pheno_df(variant)
  fixture_cached(paste0("pheno_", variant, "_", format), function(dir) {
    path <- file.path(dir, paste0("pheno.", format))
    switch(format,
      csv = utils::write.csv(df, path, row.names = FALSE),
      tsv = ,
      txt = utils::write.table(df, path, sep = "\t", row.names = FALSE, quote = FALSE),
      xlsx = {
        skip_if_not_installed("writexl")
        writexl::write_xlsx(df, path)
      },
      rds = saveRDS(df, path),
      Rdata = { pheno <- df; save(pheno, file = path) })
    path
  })
}

# A bad-sample list: one basename per line, no header.
fixture_exclude <- function() {
  fixture_cached("exclude", function(dir) {
    path <- file.path(dir, "exclude.txt")
    writeLines(fixture_basenames_450k()[2], path)
    path
  })
}

# ---- signal-intensity matrix (GEO style) ------------------------------------

# Raw (preprocessRaw) M/U intensities and detection p-values of the 6 x 450k
# samples as a tab-separated file with GEO-style headers:
#   ID_REF  <s> Unmethylated signal  <s> Methylated signal  <s> Detection Pval ...
fixture_sigmatrix <- function() {
  fixture_cached("sigmatrix", function(dir) {
    RGsetEx <- NULL
    utils::data("RGsetEx", package = "minfiData", envir = environment())
    mset <- minfi::preprocessRaw(RGsetEx)
    detP <- minfi::detectionP(RGsetEx)
    U <- minfi::getUnmeth(mset)
    M <- minfi::getMeth(mset)
    ids <- colnames(mset)
    out <- data.frame(ID_REF = rownames(mset), stringsAsFactors = FALSE)
    for (s in ids) {
      out[[paste(s, "Unmethylated signal")]] <- U[, s]
      out[[paste(s, "Methylated signal")]] <- M[, s]
      out[[paste(s, "Detection Pval")]] <- detP[rownames(mset), s]
    }
    path <- file.path(dir, "signal_matrix.tsv")
    utils::write.table(out, path, sep = "\t", row.names = FALSE, quote = FALSE)
    path
  })
}

# ---- beta matrix -------------------------------------------------------------

# Noob betas of the 6 x 450k samples with a few NAs injected. Returns a list
# with `csv` and `rds` paths and the matrix itself.
fixture_beta <- function() {
  fixture_cached("beta", function(dir) {
    RGsetEx <- NULL
    utils::data("RGsetEx", package = "minfiData", envir = environment())
    beta <- minfi::getBeta(minfi::preprocessNoob(RGsetEx, dyeMethod = "single"))
    set.seed(20260914)
    na_idx <- cbind(sample(nrow(beta), 200), sample(ncol(beta), 200, replace = TRUE))
    beta[na_idx] <- NA
    csv <- file.path(dir, "beta.csv")
    rds <- file.path(dir, "beta.rds")
    utils::write.csv(beta, csv)
    saveRDS(beta, rds)
    list(csv = csv, rds = rds, beta = beta, na_idx = na_idx)
  })
}

# ---- a run with qc_read() done ----------------------------------------------

# The 450k fixture read once per session (batch = "dir": one batch per chip);
# tests copy it with copy_run_450k() so the slow qc_read() runs only once.
fixture_run_450k <- function() {
  fixture_cached("run_450k", function(dir) {
    run_dir <- file.path(dir, "run")
    cfg <- suppressMessages(qc_config(input = fixture_idat_450k(), run_dir = run_dir,
                                      array = "450k", batch = "dir", save_rs = TRUE))
    qc_init(cfg)
    qc_read(run_dir)
    run_dir
  })
}

# A fresh copy of fixture_run_450k() in a temp dir. Extra arguments re-create
# the config (qc_init(force = TRUE)); batch must stay "dir".
copy_run_450k <- function(..., env = parent.frame()) {
  skip_if_no_fixture_data()
  src <- fixture_run_450k()
  dst <- file.path(withr::local_tempdir(.local_envir = env), "run")
  dir.create(dst)
  file.copy(list.files(src, full.names = TRUE), dst, recursive = TRUE)
  args <- list(...)
  if (length(args) > 0) {
    cfg <- suppressMessages(do.call(qc_config, c(
      list(input = fixture_idat_450k(), run_dir = dst, array = "450k", batch = "dir", save_rs = TRUE),
      args)))
    qc_init(cfg, force = TRUE)
  } else {
    cfg <- read_config(dst); cfg$run_dir <- dst; write_config(cfg, file.path(dst, "config.yml"))
  }
  dst
}

# A merged run of the 450k fixture (batch = "dir"), built once per session.
# failed_sample_threshold = 0.5: with the default 0.1 and only 6 samples every
# probe failing in a single sample is removed and nothing is left to impute.
fixture_merged_450k <- function() {
  fixture_cached("run_450k_merged", function(dir) {
    skip_if_not_installed("ChAMP"); skip_if_not_installed("impute")
    run_dir <- copy_run_450k(failed_sample_threshold = 0.5, env = teardown_env())
    qc_filter_samples(run_dir); qc_normalise(run_dir); qc_merge(run_dir)
    run_dir
  })
}

