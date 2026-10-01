# prepDNAm 0.0.0.9000

* Package skeleton. Successor to `eutopsQC`; carries over the probe lists
  (`data/`), `rhoEstimator()` and `subset_versionshared_CpGs()`.
* `normalise_bmiq()` replaces `champ.norm()` / `champ.normv2()` /
  `champ.normm()`: one array-agnostic BMIQ wrapper that takes the probe design
  from the array manifest (`minfi::getProbeType()`), so 450k, EPIC, EPIC v2
  and mouse arrays all use the same code path. Results are identical to
  `ChAMP::champ.norm()` for 450k/EPIC. Fixes a reproducibility issue where
  BMIQ output depended on the parallel backend's RNG kind; the caller's RNG
  state is now preserved. Missing values are tolerated and returned as `NA`;
  degenerate fits warn with the sample name.
* Internal array lookup table (`R/sysdata.rda`, built by `data-raw/arrays.R`)
  with `array_info()` / `supported_arrays()`: annotation packages, probe
  counts, samples per chip, probe-removal lists and feature availability for
  450k, EPIC, EPIC v2 and mouse. Replaces the hard-coded constants and
  `array == "mouse"` / `grepl("_BC11")` branches in eutopsQC (the latter
  misidentified EPIC v2 as mouse).
* `qc_config()`: one validated configuration object for a run (paths, array,
  input type inferred from `input`, pheno/exclude, batching mode, thresholds,
  feature flags coerced per array/input type, cores, seed).
* `discover_samples()`: builds the sample sheet from IDATs (recursive,
  `.idat.gz`, `_Red` check), GEO-style signal-matrix headers or beta-matrix
  column names; exact-basename join with `pheno` (missing on disk = error,
  extra on disk = dropped with a message); `exclude` from a vector, file or
  pheno column; `batch` = `"none"` / `"dir"` / `"chip"` / pheno column.
  Replaces `by.dir`, `find.files` and the `readline()` prompt.
* `qc_init()`: first stage; creates the run directory with `config.yml`,
  `samples.csv`, `run_info.json` (versions, platform, seed, checksums,
  per-stage status/timings) and `log.txt`. Refuses a run directory initialised
  with a different configuration unless `force = TRUE`.
* Stage infrastructure: every stage is idempotent (skips if its outputs
  exist), logs to `log.txt` (console too unless `options(prepDNAm.quiet = TRUE)`),
  sets the seed, and on error writes `FAILED_<stage>.txt` with the message and
  traceback before re-throwing. `read_config()`, `read_samples()`,
  `batch_dir()` exported for downstream use.
* `qc_read()`: per-batch read stage for all three input types, writing one
  bundle (`mset.rds`, `detP.rds`, `qc.rds` with `mMed`/`uMed`/`snp_outlier`,
  `rho.rds`, `sex_raw.rds`, `ctrls.rds`, `densities.rds`, optional `rs.rds`
  and `shiny.rds`; `beta_input.rds` for beta input). Checks the IDATs
  against the configured array; signal matrices are read with
  `data.table::fread()` and restricted to manifest CpG probes; M-values are
  rejected. Excluded samples are read and measured so the report can show
  them. `save_shiny` option (default off) keeps the shinyMethyl summary.
* `qc_filter_samples()`: per-batch `samples_removed.csv` with `reason`
  (`excluded` / `low_intensity` / `failed_probes`), `value` and `threshold`.
* `qc_normalise()`: BMIQ on the retained samples -> `beta_norm.rds`; adds
  `bmiq` curves to `densities.rds`. For beta input honours `normalise`
  (`TRUE` / `FALSE` / `"auto"`, the latter via a type II peak-shift heuristic,
  decision logged).
* `qc_merge()`: inner-joins the normalised batches, masks detP-failed
  values, removes probes (array lists, pooled detP failure fraction, NA
  fraction; all recorded in `probes_removed.csv`), imputes with
  `impute.knn` on the pooled matrix, calls sex on all samples (minfi's
  cutoff rule), and writes `beta_merged.rds`, `pheno_postQC.rds`,
  `sample_qc.csv`, `imputed.csv` (+ `imputed_mask.csv`), `rs_merged.rds` and
  `beta_merged_compatible.rds` (EPIC v2). Unlike eutopsQC, where the NA
  threshold and kNN imputation ran per plate, the result is independent of
  batching (tested: `batch = "dir"` and `"none"` give identical output).
* `run_prepDNAm()`: runs all stages from a config or resumes from a run
  directory; dependencies checked up front. `run_summary()` /
  `print()` give sample and probe counts, stage status and timings.
* `inst/scripts/run_prepDNAm.R`: command-line entry point taking a YAML file
  with the `qc_config()` arguments (`--config file.yml [--force] [--quiet]`).
* `qc_report()`: renders the HTML report from a run directory into
  `run_dir/report/` (rmarkdown site; pages in `inst/rmd/`). All page logic
  lives in exported helpers usable in custom reports: `qc_data()`,
  `qc_table_batches()` / `qc_table_samples()` / `qc_table_probes()` /
  `qc_thresholds()`, `qc_plot_intensities()` / `qc_plot_failed_probes()` /
  `qc_plot_controls()` / `qc_plot_densities()` / `qc_plot_rho()` /
  `qc_plot_sex()` / `qc_table_sex()`, `qc_snp_agreement()`,
  `qc_estimate_ic()`, `qc_estimate_age()`, `qc_index_cpgs()`, `qc_dimred()`
  / `qc_plot_dimred()` / `qc_plot_pc_association()`. New in the report
  compared with eutopsQC: an overview dashboard, a searchable per-sample table
  with all metrics and removal reasons, failed-probe fractions, a sex check
  against reported sex, genotype-based duplicate / swap detection
  (`ewastools`), and a PC-vs-variable association screen; the AHRR smoking CpG
  (cg05575921) replaces the misattributed ELOVL2 CpG, which is now shown as
  an age marker.
* Report revisions: every page is a self-contained HTML file by default
  (`qc_report(self_contained = FALSE)` for one shared `site_libs/`); the
  overview has a run details table (name, date, time, array, input type,
  version) and plain stat tiles; shorter navbar labels that fit on one line;
  beta densities are interactive, one colour per sample, removed samples
  dashed, sample on hover (`qc_plot_densities(n_points = )` thins the
  curves); page 6 shows the full sex table, and the immune cell, age and
  marker CpG plots are exported functions (`qc_plot_ic()`, `qc_plot_age()`,
  `qc_plot_index_cpgs()`); `qc_report()` takes `sex_col` / `donor_col`;
  PCA/UMAP plots are square and laid out in a grid; the PC association
  heatmap uses viridis. Fixed: plots printed in loops (control probes,
  densities) rendered blank because plotly.js was not included.
* `qc_plot_snp_pairs()` / `qc_snp_pairs_pdf()`: SNP probe betas of each
  donor's samples plotted against each other, as ggplots or a PDF with one
  page per donor (kept out of the HTML report, which would grow with the
  number of donors).
* Donor identifiers: the sample sheet column `donor_id` is the convention;
  common alternatives (`person`, `patient`, ...) are still recognised, with a
  message, and `donor_col` names any other column.
* `save_rs` now defaults to saving the SNP probe betas whenever the array and
  input type have them.
* `export_beta_csv()`: writes the merged betas and sample metadata as CSV in
  the DNA Methylation Array Data Standard V-2410, the format read by
  `biolearn` (`cpgSite` + one column per sample, 3 decimals, `NaN` for
  missing, at most 1000 samples per file; metadata with `SampleID`, `Sex`
  0/1, `Age`, `PredictedSex`). Probe IDs as in `beta_merged.rds` (EPIC v2
  suffixes kept); `digits = NULL` / `qc_config(export_digits = "full")`
  writes full precision instead of 3 decimals; `sample_id` picks the ID column (e.g. GEO accessions);
  `compatible = TRUE` exports the EPIC v1/v2 compatible matrix of an EPIC v2
  run instead. `qc_config(export_csv = TRUE)` runs the export at the end of
  `run_prepDNAm()` (for EPIC v2 both the full and the compatible matrix); it
  can be switched on for an existing run without re-running earlier stages.
* Sample sheet conventions (`?sample_sheet`, `vignette("sample-sheet")`):
  standard columns `sex` (F/M), `age` (years, numeric) and `donor_id`,
  validated when the run is initialised (unrecognised sex values, numeric sex
  codes, non-numeric or negative ages are an error listing the samples; ages
  above 120 warn; text ages are converted). Common alternative names
  (`gender`, `person`, `age_at_sampling`, ...) are recognised with a message.
* Documentation: README with installation and quick start;
  `vignette("prepDNAm")` (installation, configuration, running and resuming,
  cluster use, outputs, report, export) and `vignette("sample-sheet")`.
* Fixed: EPIC v2 runs failed at the signal-to-noise step ("`probe_type` must
  have one entry per row of `x`"). Probe types were read with
  `minfi::getProbeType()`, which uses the annotation package; the EPIC v2
  annotation (20a1) lacks 6,915 of the manifest's probes. Probe types now come
  from the manifest, matched by probe name (signal-to-noise and BMIQ).
* Fixed: the report failed on the identity page when `EpiDISH` or
  `WIDclocks` was installed ("object 'centEpiFibIC.m' not found"). The
  EpiDISH reference is now loaded with `data()`, and the WID general clock is
  computed from the WIDclocks coefficients directly (identical results;
  `WIDclocks::WID_clocks()` only works with WIDclocks and EpiDISH attached).
* Report: hierarchical EpiDISH immune subtypes (`qc_estimate_ic(subtypes =
  TRUE)`, `qc_plot_ic(subtypes = TRUE)`: 12 immune cell types with
  `cent12CT.m` within the immune compartment); a sex check summary (match /
  mismatch count and reported-vs-predicted table); static control probe
  plots (smaller page); PCA/UMAP plots with a fixed-size square plotting
  area (fixed legend space, legend labels truncated to 25 characters, no
  legend for variables with more than `max_legend` values).
* `qc_dimred()` / `qc_plot_pc_association()`: variables are split into
  technical (batch, Sentrix ID and position, which are derived from IDAT
  basenames if missing, plate, well, ...; `technical =` adds others) and
  biological, shown as two blocks with rows clustered within each; per-sample
  QC metrics are no longer included; smaller rows.
* `qc_estimates()`: new stage after the merge (`qc_config(estimates = TRUE)`,
  the default) that adds EpiDISH fractions (`Epi`, `Fib`, `IC` and 12 immune
  subtypes), `wid_general_age` and the marker CpG betas to
  `pheno_postQC.rds` and `estimates.csv`; the CSV export metadata includes
  them. Can be run on existing runs without repeating the merge.
