#' Report plots
#'
#' `ggplot2` objects behind the report pages; usable directly for custom
#' reports (wrap in `plotly::ggplotly()` for hover text, which is carried in
#' the `text` aesthetic). Each returns `NULL` when the run lacks the data it
#' needs (e.g. control probes for signal-matrix input).
#'
#' * `qc_plot_intensities()`: median methylated vs unmethylated log2 intensity
#'   per sample, threshold lines, removed samples marked.
#' * `qc_plot_failed_probes()`: fraction of failed probes per sample.
#' * `qc_plot_controls()`: intensities of one Illumina control probe type
#'   (`qc_control_types()` lists those available); `NORM_*` types as boxplots.
#' * `qc_plot_densities()`: per-sample beta densities at one `stage`
#'   (`"raw"`, `"noob"`, `"bmiq"`; `qc_density_stages()` lists those
#'   available), one colour per sample, removed samples dashed; hover text
#'   gives the sample.
#' * `qc_plot_rho()`: signal-to-noise estimate per sample.
#' * `qc_plot_sex()`: chrX vs chrY median log2 copy number, predicted sex and
#'   mismatches against a reported sex column (`qc_table_sex()`).
#'
#' @param x Run directory or [qc_data()] object.
#' @param type Control probe type, one of `qc_control_types(x)`.
#' @param stage Density stage, one of `qc_density_stages(x)`.
#' @param batch Batch label(s) to show; `NULL` for all (facetted).
#' @param n_points Number of points per density curve (at most the 512
#'   stored); fewer keeps interactive plots of large runs responsive. `NULL`
#'   for all.
#' @param sex_col Sample sheet column with reported sex (`F`/`M`, `female`/
#'   `male`, case-insensitive); `NULL` to auto-detect a column named `sex`.
#' @return A `ggplot` (or a character vector / data frame for the helper
#'   listings), or `NULL` if not available.
#' @examples
#' \dontrun{
#' qc_plot_intensities("runs/study1")
#' plotly::ggplotly(qc_plot_controls("runs/study1", "BISULFITE CONVERSION I"))
#' }
#' @name qc_plots
NULL

# geom_point that accepts a `text` aesthetic (hover text for ggplotly)
# without ggplot2's "unknown aesthetic" warning.
geom_point_tip <- function(...) suppressWarnings(ggplot2::geom_point(...))

theme_qc <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "bottom", strip.text = ggplot2::element_text(face = "bold"))
}

scale_removed <- function() {
  list(
    ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4), name = "removed",
                                labels = c(`FALSE` = "retained", `TRUE` = "removed")),
    ggplot2::scale_alpha_manual(values = c(`FALSE` = 0.8, `TRUE` = 1), guide = "none"))
}

# Samples with metrics, `removed` and hover text.
sample_frame <- function(rd) {
  s <- qc_table_samples(rd)
  s$removed <- factor(isTRUE_vec(s$removed), levels = c("FALSE", "TRUE"))
  s$hover <- paste0(s$basename, "\nbatch: ", s$batch,
                    ifelse(nzchar(s$reasons), paste0("\nremoved: ", s$reasons), ""))
  s
}

isTRUE_vec <- function(x) !is.na(x) & as.logical(x)

#' @rdname qc_plots
#' @export
qc_plot_intensities <- function(x) {
  rd <- qc_data(x)
  s <- sample_frame(rd)
  if (!"mMed" %in% names(s) || all(is.na(s$mMed))) return(NULL)
  thr <- rd$config$intensity_threshold
  lim <- range(c(s$mMed, s$uMed, thr), na.rm = TRUE) + c(-0.5, 0.5)
  ggplot2::ggplot(s, ggplot2::aes(x = .data$mMed, y = .data$uMed, colour = .data$batch,
                                  shape = .data$removed, alpha = .data$removed, text = .data$hover)) +
    ggplot2::geom_hline(yintercept = thr, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_vline(xintercept = thr, linetype = "dashed", colour = "grey50") +
    geom_point_tip(size = 2) +
    scale_removed() +
    ggplot2::coord_equal(xlim = lim, ylim = lim) +
    ggplot2::labs(x = "median methylated log2 intensity", y = "median unmethylated log2 intensity",
                  colour = "batch") +
    theme_qc()
}

#' @rdname qc_plots
#' @export
qc_plot_failed_probes <- function(x) {
  rd <- qc_data(x)
  s <- sample_frame(rd)
  if (!"frac_failed_probes" %in% names(s) || all(is.na(s$frac_failed_probes))) return(NULL)
  s$index <- stats::ave(seq_len(nrow(s)), s$batch, FUN = seq_along)
  thr <- rd$config$failed_probe_threshold
  ggplot2::ggplot(s, ggplot2::aes(x = .data$index, y = .data$frac_failed_probes, colour = .data$batch,
                                  shape = .data$removed, alpha = .data$removed, text = .data$hover)) +
    ggplot2::geom_hline(yintercept = thr, linetype = "dashed", colour = "grey50") +
    geom_point_tip(size = 2) +
    scale_removed() +
    ggplot2::scale_y_log10() +
    ggplot2::facet_wrap(~batch, scales = "free_x") +
    ggplot2::labs(x = "sample", y = sprintf("fraction of probes with detP > %s", rd$config$detp_threshold),
                  colour = "batch") +
    theme_qc() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
}

# Illumina control probe types in display order.
control_type_order <- c("STAINING", "EXTENSION", "TARGET REMOVAL", "BISULFITE CONVERSION I",
                        "BISULFITE CONVERSION II", "HYBRIDIZATION", "SPECIFICITY I", "SPECIFICITY II",
                        "NON-POLYMORPHIC", "NEGATIVE", "NORM_A", "NORM_C", "NORM_G", "NORM_T", "NORM_ACGT",
                        "RESTORATION")

#' @rdname qc_plots
#' @export
qc_control_types <- function(x) {
  rd <- qc_data(x)
  if (!rd$has[["ctrls"]]) return(character())
  types <- unique(unlist(lapply(batch_items(rd, "ctrls.rds"), function(c) c$info$Type)))
  c(intersect(control_type_order, types), setdiff(types, control_type_order))
}

# Long data frame of control probe intensities for one type.
control_frame <- function(rd, type) {
  parts <- lapply(rd$batches, function(b) {
    c <- readRDS(file.path(batch_dir(rd$run_dir, b), "ctrls.rds"))
    keep <- c$info$Type == type
    if (!any(keep)) return(NULL)
    info <- c$info[keep, , drop = FALSE]
    long <- function(m, fluor) {
      m <- m[info$Address, , drop = FALSE]
      data.frame(address = rep(info$Address, ncol(m)),
                 probe = rep(info$ExtendedType, ncol(m)),
                 basename = rep(colnames(m), each = nrow(m)),
                 fluorophore = fluor, intensity = as.vector(m), batch = b, stringsAsFactors = FALSE)
    }
    rbind(long(c$green, "Green"), long(c$red, "Red"))
  })
  do.call(rbind, parts)
}

#' @rdname qc_plots
#' @export
qc_plot_controls <- function(x, type) {
  rd <- qc_data(x)
  if (!rd$has[["ctrls"]]) return(NULL)
  d <- control_frame(rd, type)
  if (is.null(d) || nrow(d) == 0) return(NULL)
  s <- sample_frame(rd)
  d$removed <- s$removed[match(d$basename, s$basename)]
  d$basename <- factor(d$basename, levels = s$basename)
  d$hover <- paste0(d$basename, "\n", d$probe, "\n", round(d$intensity))
  p <- if (startsWith(type, "NORM")) {
    ggplot2::ggplot(d, ggplot2::aes(x = .data$basename, y = .data$intensity, fill = .data$batch)) +
      ggplot2::geom_boxplot(outlier.size = 0.5, linewidth = 0.3) +
      ggplot2::labs(fill = "batch")
  } else {
    ggplot2::ggplot(d, ggplot2::aes(x = .data$basename, y = .data$intensity, colour = .data$probe,
                                    shape = .data$removed, text = .data$hover)) +
      geom_point_tip(alpha = 0.8, size = 1.8) +
      ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4), guide = "none") +
      ggplot2::scale_colour_viridis_d(name = NULL)
  }
  p + ggplot2::geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey40") +
    ggplot2::facet_wrap(~fluorophore) +
    ggplot2::labs(title = type, x = "sample", y = "intensity") +
    theme_qc() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
}

#' @rdname qc_plots
#' @export
qc_density_stages <- function(x) {
  rd <- qc_data(x)
  if (!rd$has[["densities"]]) return(character())
  d <- batch_items(rd, "densities.rds")
  stages <- c("raw", "noob", "bmiq")
  stages[vapply(stages, function(s) all(vapply(d, function(z) !is.null(z[[s]]), logical(1))), logical(1))]
}

#' @rdname qc_plots
#' @export
qc_plot_densities <- function(x, stage = "bmiq", batch = NULL, n_points = NULL) {
  rd <- qc_data(x)
  if (!stage %in% qc_density_stages(rd)) return(NULL)
  batches <- if (is.null(batch)) rd$batches else batch
  s <- sample_frame(rd)
  parts <- lapply(batches, function(b) {
    d <- readRDS(file.path(batch_dir(rd$run_dir, b), "densities.rds"))
    m <- d[[stage]]
    if (is.null(m) || ncol(m) == 0) return(NULL)
    keep <- seq_len(nrow(m))
    if (!is.null(n_points) && n_points < nrow(m)) keep <- unique(round(seq(1, nrow(m), length.out = n_points)))
    m <- m[keep, , drop = FALSE]
    data.frame(x = rep(d$x[keep], ncol(m)), density = as.vector(m),
               basename = rep(colnames(m), each = nrow(m)), batch = b, stringsAsFactors = FALSE)
  })
  long <- do.call(rbind, parts)
  if (is.null(long)) return(NULL)
  i <- match(long$basename, s$basename)
  long$removed <- factor(s$removed[i], levels = c("FALSE", "TRUE"))
  long$hover <- s$hover[i]
  long$basename <- factor(long$basename, levels = unique(long$basename))
  # one hue per sample (too many for a legend: identify samples by hover)
  ggplot2::ggplot(long, ggplot2::aes(x = .data$x, y = .data$density, group = .data$basename,
                                     colour = .data$basename, linetype = .data$removed,
                                     text = .data$hover)) +
    suppressWarnings(ggplot2::geom_line(alpha = 0.7, linewidth = 0.4)) +
    ggplot2::scale_colour_hue(l = 55, guide = "none") +
    ggplot2::scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "22"), guide = "none") +
    ggplot2::facet_wrap(~batch) +
    ggplot2::labs(x = "beta", y = "density", title = paste("Beta densities:", stage)) +
    theme_qc() +
    ggplot2::theme(legend.position = "none")
}

#' @rdname qc_plots
#' @export
qc_plot_rho <- function(x) {
  rd <- qc_data(x)
  s <- sample_frame(rd)
  if (!"rho" %in% names(s) || all(is.na(s$rho))) return(NULL)
  s$index <- stats::ave(seq_len(nrow(s)), s$batch, FUN = seq_along)
  ggplot2::ggplot(s, ggplot2::aes(x = .data$index, y = .data$rho, colour = .data$batch,
                                  shape = .data$removed, alpha = .data$removed, text = .data$hover)) +
    geom_point_tip(size = 2) +
    scale_removed() +
    ggplot2::facet_wrap(~batch, scales = "free_x") +
    ggplot2::labs(x = "sample", y = "signal-to-noise (rho)", colour = "batch") +
    theme_qc() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
}

#' @rdname qc_plots
#' @export
qc_table_sex <- function(x, sex_col = NULL) {
  rd <- qc_data(x)
  s <- sample_frame(rd)
  if (!"predicted_sex" %in% names(s) || all(is.na(s$predicted_sex))) return(NULL)
  col <- find_sex_col(s, sex_col)
  s$reported_sex <- if (is.null(col)) NA_character_ else normalise_sex(s[[col]])
  s$sex_mismatch <- !is.na(s$reported_sex) & !is.na(s$predicted_sex) & s$reported_sex != s$predicted_sex
  s[, c("basename", "batch", "xMed", "yMed", "predicted_sex", "reported_sex", "sex_mismatch", "removed")]
}

#' @rdname qc_plots
#' @export
qc_plot_sex <- function(x, sex_col = NULL) {
  t <- qc_table_sex(x, sex_col)
  if (is.null(t)) return(NULL)
  t$hover <- paste0(t$basename, "\npredicted: ", t$predicted_sex,
                    ifelse(is.na(t$reported_sex), "", paste0("\nreported: ", t$reported_sex)))
  t$status <- ifelse(t$sex_mismatch, "mismatch", ifelse(is.na(t$reported_sex), "no reported sex", "match"))
  lim <- range(c(t$xMed, t$yMed), na.rm = TRUE)
  ggplot2::ggplot(t, ggplot2::aes(x = .data$xMed, y = .data$yMed, colour = .data$predicted_sex,
                                  shape = .data$status, text = .data$hover)) +
    ggplot2::geom_abline(intercept = -2, slope = 1, linetype = "dashed", colour = "grey50") +
    geom_point_tip(size = 2.2, alpha = 0.85) +
    ggplot2::scale_colour_manual(values = c(F = "#e7298a", M = "#1b9e77"), name = "predicted sex",
                                 na.value = "grey60") +
    ggplot2::scale_shape_manual(values = c(match = 16, mismatch = 17, `no reported sex` = 1), name = NULL) +
    ggplot2::coord_equal(xlim = lim + c(-0.3, 0.3), ylim = lim + c(-0.3, 0.3)) +
    ggplot2::labs(x = "chrX median log2 copy number", y = "chrY median log2 copy number") +
    theme_qc()
}
