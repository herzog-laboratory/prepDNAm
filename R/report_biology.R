#' Sample identity and biology checks
#'
#' Helpers behind report page 5, all working on the retained samples of a
#' finished run. Each returns `NULL` (with a message) when the data or an
#' optional package is not available.
#'
#' * `qc_snp_agreement()`: genotype calls from the SNP (`rs`) probes via
#'   `ewastools`. Samples are grouped into inferred donors
#'   (`ewastools::enumerate_sample_donors()`); if the sample sheet has a donor
#'   column (`donor_col`; by convention `donor_id`), samples whose genotypes
#'   disagree with their reported donor, or agree with another donor's, are
#'   listed (`ewastools::check_snp_agreement()`). Without `donor_id`, common
#'   alternatives (`donor`, `person`, `person_id`, `individual`, `subject`,
#'   `patient`, `participant`) are used with a message.
#' * `qc_estimate_ic()`: epithelial / fibroblast / immune cell fractions
#'   (`EpiDISH`, robust partial correlations with `centEpiFibIC.m`); with
#'   `subtypes = TRUE` the immune fraction is split into 12 immune cell types
#'   (hierarchical EpiDISH with `cent12CT.m`).
#' * `qc_estimate_age()`: WID general clock age (`WIDclocks`).
#' * `qc_index_cpgs()`: betas of two well-characterised marker CpGs:
#'   `cg05575921` (AHRR, hypomethylated in smokers) and `cg16867657` (ELOVL2,
#'   increases with age).
#'
#' * `qc_plot_ic()`, `qc_plot_age()`, `qc_plot_index_cpgs()`: plots of the
#'   estimates above (`NULL` when the estimate is not available). Age is
#'   plotted against chronological age if the sample sheet has an age column
#'   (`age_col`, auto-detected among `age`, `age_at_sampling`, `chron_age`).
#' * `qc_plot_snp_pairs()`: for every donor with more than one sample, SNP
#'   probe betas of each pair of its samples against each other (one panel per
#'   pair, Pearson r in the panel title). Samples of the same donor lie on the
#'   diagonal in three clusters (genotypes AA/AB/BB); a scattered panel points
#'   to a swapped or mislabelled sample. Donors come from `donor_col` if the
#'   sample sheet has one, otherwise from the genotype-inferred donors.
#'   `qc_snp_pairs_pdf()` writes these plots to a PDF, one page per donor
#'   (not part of the HTML report, which would grow with the number of
#'   donors).
#'
#' Indices trained on EPIC v1 (`EpiDISH`, `WIDclocks`, marker CpGs) use the
#' EPIC v1/v2 compatible matrix for EPIC v2 runs. Human arrays only.
#'
#' @param x Run directory or [qc_data()] object.
#' @param subtypes `qc_estimate_ic()` / `qc_plot_ic()`: also split the immune
#'   compartment into 12 immune cell types with hierarchical EpiDISH
#'   (`EpiDISH::hepidish()`, `centEpiFibIC.m` then `cent12CT.m` within the
#'   immune compartment). The subtype reference is derived from blood, so
#'   treat subtype fractions in solid tissues as indicative.
#' @param donor_col Sample sheet column with donor identifiers; `NULL` to
#'   auto-detect.
#' @param age_col Sample sheet column with chronological age; `NULL` to
#'   auto-detect.
#' @param donors Donors to plot (values of the donor column or inferred donor
#'   numbers); `NULL` for all with more than one sample.
#' @param max_pairs Donors with more sample pairs than this are shown as each
#'   sample against the donor's first sample instead of all pairs.
#' @param file PDF to write; defaults to `snp_pairs.pdf` in the run
#'   directory.
#' @return `qc_snp_agreement()`: a list with `donors` (per sample:
#'   `inferred_donor`, `n_samples_of_donor`, reported donor if known),
#'   `conflicts` (data frame, possibly empty; `NULL` without a donor column)
#'   and `donor_col`. `qc_estimate_*()` / `qc_index_cpgs()`: a data frame with
#'   `basename` and the estimates. `qc_plot_*()`: a `ggplot`, or `NULL`;
#'   `qc_plot_snp_pairs()` a named list of `ggplot`s, one per donor;
#'   `qc_snp_pairs_pdf()` the path of the PDF, invisibly (`NULL` if there was
#'   nothing to plot).
#' @examples
#' \dontrun{
#' qc_snp_agreement("runs/study1")$donors
#' qc_plot_age("runs/study1")
#' qc_plot_snp_pairs("runs/study1", donor_col = "patient_id")[[1]]
#' qc_snp_pairs_pdf("runs/study1")
#' }
#' @name qc_biology
NULL

#' @rdname qc_biology
#' @export
qc_snp_agreement <- function(x, donor_col = NULL) {
  rd <- qc_data(x)
  if (!rd$has[["rs"]]) {
    rlang::inform("No rs_merged.rds: SNP probes were not saved (`save_rs = TRUE`) or are not available.")
    return(NULL)
  }
  if (!requireNamespace("ewastools", quietly = TRUE)) {
    rlang::inform("Package ewastools is not installed; genotype checks skipped.")
    return(NULL)
  }
  rs <- readRDS(run_file(rd$run_dir, "rs_merged.rds"))
  if (ncol(rs) < 2) return(NULL)
  genotypes <- ewastools::call_genotypes(rs, learn = FALSE)
  inferred <- ewastools::enumerate_sample_donors(genotypes)
  donors <- data.frame(basename = colnames(rs), inferred_donor = as.integer(inferred),
                       stringsAsFactors = FALSE)
  donors$n_samples_of_donor <- as.integer(table(donors$inferred_donor)[as.character(donors$inferred_donor)])
  s <- qc_table_samples(rd)
  col <- find_donor_col(s, donor_col)
  conflicts <- NULL
  if (!is.null(col)) {
    donors$reported_donor <- as.character(s[[col]][match(donors$basename, s$basename)])
    ok <- !is.na(donors$reported_donor)
    res <- ewastools::check_snp_agreement(genotypes, donor_ids = donors$reported_donor,
                                          sample_ids = donors$basename)
    conflicts <- if (is.null(res) || length(res) == 0) {
      data.frame(donor1 = character(), sample1 = character(), donor2 = character(),
                 sample2 = character(), agreement = numeric(), stringsAsFactors = FALSE)
    } else {
      as.data.frame(do.call(rbind, lapply(res, as.data.frame)), stringsAsFactors = FALSE)
    }
    rownames(conflicts) <- NULL
  }
  list(donors = donors, conflicts = conflicts, donor_col = col)
}

human_only <- function(rd, what) {
  if (rd$array$species != "human" || !rd$array$has_clocks) {
    rlang::inform(sprintf("%s is only available for human arrays.", what))
    return(FALSE)
  }
  if (!rd$has[["beta_merged"]]) {
    rlang::inform("No beta_merged.rds; run qc_merge() first.")
    return(FALSE)
  }
  TRUE
}

#' @rdname qc_biology
#' @export
qc_estimate_ic <- function(x, subtypes = FALSE) {
  rd <- qc_data(x)
  check_flag(subtypes, "subtypes")
  if (!human_only(rd, "Immune cell estimation")) return(NULL)
  if (!requireNamespace("EpiDISH", quietly = TRUE)) {
    rlang::inform("Package EpiDISH is not installed; immune cell estimation skipped.")
    return(NULL)
  }
  beta <- report_beta(rd, compatible = TRUE)
  est <- if (subtypes) {
    EpiDISH::hepidish(beta.m = beta, ref1.m = epidish_ref("centEpiFibIC.m"),
                      ref2.m = epidish_ref("cent12CT.m"), h.CT.idx = 3, method = "RPC")
  } else {
    EpiDISH::epidish(beta.m = beta, ref.m = epidish_ref("centEpiFibIC.m"), method = "RPC")$estF
  }
  data.frame(basename = rownames(est), as.data.frame(est), row.names = NULL,
             stringsAsFactors = FALSE, check.names = FALSE)
}

# An EpiDISH reference matrix: datasets of the package (data()), not objects
# in its namespace.
epidish_ref <- function(name) {
  e <- new.env()
  utils::data(list = name, package = "EpiDISH", envir = e)
  e[[name]]
}

#' @rdname qc_biology
#' @export
qc_estimate_age <- function(x) {
  rd <- qc_data(x)
  if (!human_only(rd, "Epigenetic age")) return(NULL)
  if (!requireNamespace("WIDclocks", quietly = TRUE)) {
    rlang::inform("Package WIDclocks is not installed; epigenetic age skipped.")
    return(NULL)
  }
  beta <- report_beta(rd, compatible = TRUE)
  # The WID general clock is linear: intercept + sum(beta * weight). Computed
  # here from the WIDclocks coefficients, as WIDclocks::WID_clocks() does;
  # that function finds its data and EpiDISH only when both are attached.
  e <- new.env()
  utils::data(list = "WID_general_clock_coef", package = "WIDclocks", envir = e)
  w <- e$WID_general_clock_coef
  intercept <- w[1]
  w <- w[-1]
  cpgs <- intersect(names(w), rownames(beta))
  if (length(cpgs) < length(w)) {
    rlang::inform(sprintf("WID general clock: %d of %d CpGs in the beta matrix.", length(cpgs), length(w)))
  }
  age <- intercept + colSums(beta[cpgs, , drop = FALSE] * w[cpgs])
  data.frame(basename = colnames(beta), wid_general_age = as.numeric(age), stringsAsFactors = FALSE)
}

#' @rdname qc_biology
#' @export
qc_index_cpgs <- function(x) {
  rd <- qc_data(x)
  if (!human_only(rd, "Marker CpGs")) return(NULL)
  beta <- report_beta(rd, compatible = TRUE)
  cpgs <- c(cg05575921 = "AHRR (smoking)", cg16867657 = "ELOVL2 (age)")
  present <- names(cpgs)[names(cpgs) %in% rownames(beta)]
  if (length(present) == 0) {
    rlang::inform("Marker CpGs cg05575921 / cg16867657 are not in the beta matrix.")
    return(NULL)
  }
  out <- data.frame(basename = colnames(beta), stringsAsFactors = FALSE)
  for (cg in present) out[[cg]] <- as.numeric(beta[cg, ])
  attr(out, "labels") <- cpgs[present]
  out
}

#' @rdname qc_biology
#' @export
qc_plot_ic <- function(x, subtypes = FALSE) {
  ic <- qc_estimate_ic(x, subtypes = subtypes)
  if (is.null(ic)) return(NULL)
  types <- setdiff(names(ic), "basename")
  long <- data.frame(basename = rep(ic$basename, length(types)), type = rep(types, each = nrow(ic)),
                     fraction = unlist(ic[types], use.names = FALSE), stringsAsFactors = FALSE)
  long$basename <- factor(long$basename, levels = ic$basename)
  long$type <- factor(long$type, levels = rev(types))    # Epi at the bottom of each bar
  long$hover <- paste0(long$basename, "\n", long$type, ": ", round(long$fraction, 3))
  fill <- if (subtypes) {
    # Epi / Fib in greys, the immune subtypes in colour
    immune <- setdiff(types, c("Epi", "Fib"))
    ggplot2::scale_fill_manual(
      values = c(Epi = "grey80", Fib = "grey55",
                 stats::setNames(grDevices::hcl.colors(length(immune), "Dark 3"), immune)),
      breaks = types, name = NULL)
  } else {
    ggplot2::scale_fill_brewer(palette = "Set2", breaks = types, name = NULL)
  }
  ggplot2::ggplot(long, ggplot2::aes(x = .data$basename, y = .data$fraction, fill = .data$type,
                                     text = .data$hover)) +
    suppressWarnings(ggplot2::geom_col()) +
    fill +
    ggplot2::labs(x = "sample", y = "fraction") +
    theme_qc() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
}

#' @rdname qc_biology
#' @export
qc_plot_age <- function(x, age_col = NULL) {
  rd <- qc_data(x)
  s <- qc_table_samples(rd)
  col <- find_age_col(s, age_col)
  age <- qc_estimate_age(rd)
  if (is.null(age)) return(NULL)
  age$age <- if (is.null(col)) NA_real_ else as.numeric(s[[col]][match(age$basename, s$basename)])
  if (all(is.na(age$age))) {
    age$index <- seq_len(nrow(age))
    age$hover <- paste0(age$basename, "\nWID age: ", round(age$wid_general_age, 1))
    p <- ggplot2::ggplot(age, ggplot2::aes(x = .data$index, y = .data$wid_general_age, text = .data$hover)) +
      geom_point_tip(size = 2, alpha = 0.8, colour = "grey30") +
      ggplot2::labs(x = "sample", y = "WID general age")
  } else {
    age$hover <- paste0(age$basename, "\nage: ", age$age, "\nWID age: ", round(age$wid_general_age, 1))
    lim <- range(c(age$age, age$wid_general_age), na.rm = TRUE)
    p <- ggplot2::ggplot(age, ggplot2::aes(x = .data$age, y = .data$wid_general_age, text = .data$hover)) +
      ggplot2::geom_abline(linetype = "dashed", colour = "grey60") +
      geom_point_tip(size = 2, alpha = 0.8, colour = "grey30") +
      ggplot2::coord_equal(xlim = lim, ylim = lim) +
      ggplot2::labs(x = sprintf("chronological age (%s)", col), y = "WID general age")
  }
  p + theme_qc()
}

#' @rdname qc_biology
#' @export
qc_plot_index_cpgs <- function(x) {
  mk <- qc_index_cpgs(x)
  if (is.null(mk)) return(NULL)
  labs <- attr(mk, "labels")
  long <- do.call(rbind, lapply(names(labs), function(cg) {
    data.frame(basename = mk$basename, cpg = sprintf("%s %s", cg, labs[[cg]]), beta = mk[[cg]],
               stringsAsFactors = FALSE)
  }))
  long$basename <- factor(long$basename, levels = mk$basename)
  long$hover <- paste0(long$basename, "\n", long$cpg, ": ", round(long$beta, 3))
  ggplot2::ggplot(long, ggplot2::aes(x = .data$basename, y = .data$beta, colour = .data$cpg,
                                     text = .data$hover)) +
    geom_point_tip(size = 2) +
    ggplot2::facet_wrap(~cpg) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = "sample", y = "beta") +
    theme_qc() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
                   legend.position = "none")
}

#' @rdname qc_biology
#' @export
qc_plot_snp_pairs <- function(x, donor_col = NULL, donors = NULL, max_pairs = 10) {
  rd <- qc_data(x)
  g <- qc_snp_agreement(rd, donor_col = donor_col)
  if (is.null(g)) return(NULL)
  rs <- readRDS(run_file(rd$run_dir, "rs_merged.rds"))
  d <- g$donors
  d$donor <- if (!is.null(g$donor_col)) d$reported_donor else as.character(d$inferred_donor)
  d <- d[!is.na(d$donor), , drop = FALSE]
  groups <- split(d$basename, d$donor)
  groups <- groups[lengths(groups) > 1]
  if (!is.null(donors)) groups <- groups[intersect(as.character(donors), names(groups))]
  label <- if (!is.null(g$donor_col)) g$donor_col else "inferred donor"
  lapply(stats::setNames(names(groups), names(groups)), function(donor) {
    smp <- groups[[donor]]
    pairs <- utils::combn(smp, 2, simplify = FALSE)
    if (length(pairs) > max_pairs) pairs <- lapply(smp[-1], function(b) c(smp[1], b))
    long <- do.call(rbind, lapply(pairs, function(pr) {
      a <- rs[, pr[1]]; b <- rs[, pr[2]]
      r <- suppressWarnings(stats::cor(a, b, use = "pairwise.complete.obs"))
      data.frame(probe = rownames(rs), a = a, b = b,
                 pair = sprintf("%s vs\n%s (r = %.2f)", pr[1], pr[2], r), stringsAsFactors = FALSE)
    }))
    long$pair <- factor(long$pair, levels = unique(long$pair))
    long$hover <- paste0(long$probe, "\n", round(long$a, 2), " / ", round(long$b, 2))
    ggplot2::ggplot(long, ggplot2::aes(x = .data$a, y = .data$b, text = .data$hover)) +
      ggplot2::geom_abline(linetype = "dashed", colour = "grey60") +
      geom_point_tip(size = 1.6, alpha = 0.7, colour = "grey20") +
      ggplot2::facet_wrap(~pair) +
      ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
      ggplot2::labs(x = "SNP probe beta, first sample", y = "SNP probe beta, second sample",
                    title = sprintf("%s: %s", label, donor)) +
      theme_qc()
  })
}

#' @rdname qc_biology
#' @export
qc_snp_pairs_pdf <- function(x, file = NULL, donor_col = NULL, donors = NULL, max_pairs = 10) {
  rd <- qc_data(x)
  if (is.null(file)) file <- run_file(rd$run_dir, "snp_pairs.pdf")
  plots <- qc_plot_snp_pairs(rd, donor_col = donor_col, donors = donors, max_pairs = max_pairs)
  if (length(plots) == 0) {
    rlang::inform("No donor with more than one sample: nothing to plot.")
    return(invisible(NULL))
  }
  grDevices::pdf(file, width = 10, height = 7)
  on.exit(grDevices::dev.off())
  for (p in plots) print(p)
  invisible(file)
}

#' Dimensionality reduction of the merged beta matrix
#'
#' PCA (and UMAP if `uwot` is installed) on the most variable probes of
#' `beta_merged`, plus a test of each principal component against the sample
#' sheet variables (Kruskal-Wallis for categorical, Pearson correlation for
#' numeric variables) to spot batch effects.
#'
#' Variables are the sample sheet columns with 2 to n - 1 distinct values
#' (identifiers, paths, dates and per-sample QC metrics are left out), split
#' into two groups:
#' * **technical**: `batch`, the Sentrix ID and position (parsed from IDAT
#'   basenames `<sentrix_id>_<sentrix_position>` if the sample sheet has no
#'   such columns), and columns whose name refers to processing (plate, well,
#'   slide, chip, array, position, batch, run, scan, pool, date, row, column);
#'   `technical` adds others;
#' * **biological**: all other columns.
#'
#' @param x Run directory or [qc_data()] object.
#' @param n_probes Number of most variable probes to use.
#' @param n_pcs Number of principal components to test.
#' @param seed Seed for UMAP.
#' @param technical Further sample sheet columns to treat as technical.
#' @param cluster Order the rows of each group by hierarchical clustering of
#'   their association profiles.
#' @param dr Result of `qc_dimred()`.
#' @param method `"pca"` or `"umap"`.
#' @param colour Column of `dr$coords` to colour by (see `dr$variables`).
#' @param max_legend Categorical variables with more values than this are
#'   plotted without a legend (hover shows each sample's value). Legend labels
#'   are truncated to 25 characters.
#' @return `qc_dimred()`: a list of class `prepDNAm_dimred` with `coords`
#'   (per sample: `PC1..`, `UMAP1`, `UMAP2` and the variables),
#'   `var_explained`, `pvals` (variables x PCs), `variables`, `groups` (named:
#'   variable -> `"technical"` / `"biological"`), `n_probes`, `has_umap`. `qc_plot_dimred()` / `qc_plot_pc_association()`:
#'   a `ggplot`.
#' @examples
#' \dontrun{
#' dr <- qc_dimred("runs/study1")
#' qc_plot_dimred(dr, "pca", colour = "batch")
#' qc_plot_pc_association(dr)
#' }
#' @export
qc_dimred <- function(x, n_probes = 30000, n_pcs = 10, seed = 1, technical = NULL) {
  rd <- qc_data(x)
  if (!rd$has[["beta_merged"]]) rlang::abort("No beta_merged.rds; run qc_merge() first.")
  beta <- report_beta(rd)
  if (ncol(beta) < 3) rlang::abort("Dimensionality reduction needs at least 3 samples.")
  beta <- beta[stats::complete.cases(beta), , drop = FALSE]
  sds <- matrixStats::rowSds(beta)
  top <- order(sds, decreasing = TRUE)[seq_len(min(n_probes, nrow(beta)))]
  b <- t(beta[top, , drop = FALSE])

  n_pcs <- min(n_pcs, nrow(b) - 1)
  pca <- stats::prcomp(b, center = TRUE, scale. = FALSE, rank. = n_pcs)
  coords <- data.frame(basename = rownames(b), pca$x[, seq_len(n_pcs), drop = FALSE],
                       stringsAsFactors = FALSE)
  var_explained <- stats::setNames((pca$sdev^2 / sum(pca$sdev^2))[seq_len(n_pcs)],
                                   paste0("PC", seq_len(n_pcs)))

  has_umap <- requireNamespace("uwot", quietly = TRUE)
  if (has_umap) {
    set.seed(seed)
    um <- uwot::umap(b, n_neighbors = min(15, nrow(b) - 1), min_dist = 0.1, seed = seed)
    coords$UMAP1 <- um[, 1]; coords$UMAP2 <- um[, 2]
  }

  s <- qc_table_samples(rd)
  s <- s[match(coords$basename, s$basename), , drop = FALSE]
  s <- add_sentrix_columns(s)
  qc_metrics <- setdiff(names(rd$sample_qc), c("basename", "batch"))
  vars <- setdiff(pheno_variables(s), c(qc_metrics, "removed", "reasons"))
  groups <- stats::setNames(ifelse(is_technical(vars) | vars %in% technical, "technical", "biological"), vars)
  vars <- c(vars[groups == "technical"], vars[groups == "biological"])
  groups <- groups[vars]
  for (v in vars) coords[[v]] <- s[[v]]

  pcs <- paste0("PC", seq_len(n_pcs))
  pvals <- matrix(NA_real_, nrow = length(vars), ncol = n_pcs, dimnames = list(vars, pcs))
  for (v in vars) {
    y <- s[[v]]
    for (pc in pcs) {
      pvals[v, pc] <- tryCatch({
        if (is.numeric(y)) stats::cor.test(coords[[pc]], y)$p.value
        else stats::kruskal.test(coords[[pc]], factor(y))$p.value
      }, error = function(e) NA_real_)
    }
  }
  structure(list(coords = coords, var_explained = var_explained, pvals = pvals,
                 variables = vars, groups = groups, n_probes = length(top), has_umap = has_umap),
            class = "prepDNAm_dimred")
}

# Sentrix ID / position from IDAT basenames (<id>_R01C01), unless the
# sample sheet already has them.
add_sentrix_columns <- function(s) {
  m <- regmatches(s$basename, regexec("^([0-9]+)_(R[0-9]{2}C[0-9]{2})$", s$basename))
  ok <- lengths(m) == 3
  if (!any(ok)) return(s)
  lower <- tolower(names(s))
  if (!any(lower %in% c("sentrix_id", "sentrixid", "slide", "chip"))) {
    s$sentrix_id <- ifelse(ok, vapply(m, function(z) if (length(z) == 3) z[2] else NA_character_, ""), NA)
  }
  if (!any(lower %in% c("sentrix_position", "sentrixposition", "array", "position"))) {
    s$sentrix_position <- ifelse(ok, vapply(m, function(z) if (length(z) == 3) z[3] else NA_character_, ""), NA)
  }
  s
}

# Column names that refer to sample processing rather than biology.
is_technical <- function(v) {
  grepl("batch|plate|well|sentrix|slide|chip|array|position|scan|pool|date|(^|_)(run|row|col|column)($|_)",
        tolower(v))
}

#' @rdname qc_dimred
#' @export
qc_plot_dimred <- function(dr, method = c("pca", "umap"), colour = NULL, max_legend = 20) {
  method <- match.arg(method)
  if (method == "umap" && !dr$has_umap) return(NULL)
  d <- dr$coords
  xy <- if (method == "pca") c("PC1", "PC2") else c("UMAP1", "UMAP2")
  labs <- if (method == "pca") {
    sprintf("%s (%.1f%%)", xy, 100 * dr$var_explained[1:2])
  } else xy
  d$hover <- if (is.null(colour)) d$basename else paste0(d$basename, "\n", colour, ": ", d[[colour]])
  title <- if (is.null(colour)) toupper(method) else colour
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data[[xy[1]]], y = .data[[xy[2]]], text = .data$hover))
  p <- if (is.null(colour)) {
    p + geom_point_tip(size = 2, alpha = 0.8, colour = "grey30")
  } else if (is.numeric(d[[colour]])) {
    p + geom_point_tip(ggplot2::aes(colour = .data[[colour]]), size = 2, alpha = 0.8) +
      ggplot2::scale_colour_viridis_c(name = NULL)
  } else {
    n <- length(unique(stats::na.omit(d[[colour]])))
    # many values: no legend (colours are indistinguishable anyway; hover shows the value)
    if (n > max_legend) title <- sprintf("%s (%d values; hover for each sample)", colour, n)
    p + geom_point_tip(ggplot2::aes(colour = as.character(.data[[colour]])), size = 2, alpha = 0.8) +
      ggplot2::scale_colour_discrete(name = NULL, labels = function(x) truncate_label(x, 25),
                                     guide = if (n > max_legend) "none" else "legend")
  }
  p + ggplot2::labs(x = labs[1], y = labs[2], title = title) +
    theme_qc() +
    ggplot2::theme(aspect.ratio = 1, legend.position = "right")   # square panel whatever the legend size
}

truncate_label <- function(x, n) ifelse(nchar(x) > n, paste0(substr(x, 1, n - 1), "\u2026"), x)

#' @rdname qc_dimred
#' @export
qc_plot_pc_association <- function(dr, cluster = TRUE) {
  if (length(dr$variables) == 0) return(NULL)
  m <- dr$pvals
  neglog <- -log10(pmax(m, 1e-16))
  groups <- if (is.null(dr$groups)) stats::setNames(rep("variables", nrow(m)), rownames(m)) else dr$groups
  # row order: technical then biological, clustered within each group
  order_rows <- function(vars) {
    if (!cluster || length(vars) < 3) return(vars)
    z <- neglog[vars, , drop = FALSE]
    z[is.na(z)] <- 0
    vars[stats::hclust(stats::dist(z))$order]
  }
  grp_levels <- intersect(c("technical", "biological", "variables"), unique(groups))
  rows <- unlist(lapply(grp_levels, function(g) order_rows(names(groups)[groups == g])))
  d <- data.frame(variable = rep(rownames(m), ncol(m)), pc = rep(colnames(m), each = nrow(m)),
                  p = as.vector(m), neglog = as.vector(neglog), stringsAsFactors = FALSE)
  d$group <- factor(groups[d$variable], levels = grp_levels)
  d$pc <- factor(d$pc, levels = colnames(m))
  d$variable <- factor(d$variable, levels = rev(rows))
  d$label <- ifelse(is.na(d$p), "", formatC(d$p, format = "g", digits = 2))
  # white labels on the dark end of the (reversed) viridis scale
  dark <- if (all(is.na(d$neglog))) Inf else 0.55 * max(d$neglog, na.rm = TRUE)
  ggplot2::ggplot(d, ggplot2::aes(x = .data$pc, y = .data$variable, fill = .data$neglog,
                                  text = paste0(.data$variable, " ~ ", .data$pc, "\np = ", .data$label))) +
    suppressWarnings(ggplot2::geom_tile(colour = "white")) +
    ggplot2::geom_text(ggplot2::aes(label = .data$label, colour = .data$neglog > dark), size = 2.3) +
    ggplot2::scale_fill_viridis_c(name = "-log10 p", direction = -1, na.value = "grey90") +
    ggplot2::scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "white"), guide = "none") +
    ggplot2::facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_qc() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text = ggplot2::element_text(size = 8),
                   strip.text.y = ggplot2::element_text(angle = 0, hjust = 0))
}
