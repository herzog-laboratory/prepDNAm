# Shared setup sourced by every page: locate the run directory and load it.
knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE, fig.align = "center",
                      fig.width = 9, fig.height = 5)
# options written by prepDNAm::qc_report(); params$run_dir when rendering a page by hand
opts <- if (file.exists("_report_options.yml")) yaml::read_yaml("_report_options.yml") else list()
run_dir <- if (!is.null(params$run_dir)) params$run_dir else opts$run_dir
if (is.null(run_dir)) stop("Pass `params = list(run_dir = ...)` or render via prepDNAm::qc_report().")
rd <- prepDNAm::qc_data(run_dir)
cfg <- rd$config

note <- function(...) cat("<div class=\"alert alert-secondary\" style=\"background:#f3f4f6;padding:.6em 1em;border-left:4px solid #adb5bd\">", ..., "</div>\n")
not_available <- function(what) note(sprintf("**%s** is not available for %s input on the %s array.", what, cfg$input_type, rd$array$label))
interactive_plot <- function(p, height = NULL, width = NULL) {
  if (is.null(p)) return(invisible(NULL))
  plotly::ggplotly(p, tooltip = "text", height = height, width = width) |> plotly::config(displayModeBar = FALSE)
}
# Print a widget from inside a loop in a results="asis" chunk; print() alone
# drops the widget's JavaScript dependencies, so register them explicitly.
emit_widget <- function(w) {
  if (is.null(w)) return(invisible(NULL))
  out <- knitr::knit_print(htmltools::tagList(w))
  knitr::knit_meta_add(attr(out, "knit_meta"))
  cat(out)
}
# Interactive scatter plot with a square plotting area of `size` px. Margins
# are fixed (plotly would otherwise widen them to fit the legend and shrink
# the plot); the legend gets a `legend_width` px strip on the right.
square_plot <- function(p, size = 380, legend_width = 190) {
  if (is.null(p)) return(invisible(NULL))
  pp <- ggplot2::ggplot_build(p)$layout$panel_params[[1]]
  xr <- pp$x.range; yr <- pp$y.range
  m <- list(l = 60, r = legend_width, t = 50, b = 50, autoexpand = FALSE)
  plotly::ggplotly(p, tooltip = "text", width = m$l + size + m$r, height = m$t + size + m$b) |>
    plotly::layout(
      xaxis = list(range = xr, constrain = "domain"),
      yaxis = list(range = yr, scaleanchor = "x", scaleratio = diff(xr) / diff(yr), constrain = "domain"),
      legend = list(x = 1.02, xanchor = "left", y = 1, yanchor = "top", font = list(size = 10)),
      margin = m) |>
    plotly::config(displayModeBar = FALSE)
}
# plotly widgets side by side, wrapping to the page width
plot_grid <- function(widgets) htmltools::div(style = "display:flex;flex-wrap:wrap;gap:1em", widgets)
qc_table <- function(df, round_cols = NULL, digits = 2, page = 25, ...) {
  if (is.null(df)) return(invisible(NULL))
  t <- DT::datatable(df, rownames = FALSE, filter = "top", extensions = "Buttons",
                     options = list(pageLength = page, scrollX = TRUE, dom = "Bfrtip",
                                    buttons = c("copy", "csv")), ...)
  round_cols <- intersect(round_cols, names(df))
  if (length(round_cols) > 0) t <- DT::formatRound(t, round_cols, digits)
  t
}
