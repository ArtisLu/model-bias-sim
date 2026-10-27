# Figure functions (saved by targets, used by the report).

viz_cols <- c(blue = "#2a78d6", orange = "#eb6834", ink = "#0b0b0b", ink2 = "#52514e", grid = "#e4e3df")

#' Shared minimal theme
theme_qba <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = viz_cols[["grid"]], linewidth = 0.3),
      axis.text = ggplot2::element_text(colour = viz_cols[["ink2"]]),
      axis.title = ggplot2::element_text(colour = viz_cols[["ink2"]]),
      strip.text = ggplot2::element_text(colour = viz_cols[["ink"]], face = "bold", hjust = 0),
      legend.position = "top"
    )
}

#' Distribution of replayed estimates at the calibrated effect theta*, one panel per dataset
#' @param sim_all Stacked replicates (with `dataset`); @param perf_all Stacked performance.
#' @param labels Named character vector: dataset key -> display label.
#' @return ggplot.
plot_estimates <- function(sim_all, perf_all, labels) {
  # Hazard-ratio analyses only (the survival-difference analysis has its own figures)
  po <- perf_all[scenario == "observed" & dataset %in% names(labels)]
  d <- merge(sim_all[!is.na(estimate)], po[, .(dataset, theta)], by = c("dataset", "theta"))
  lines <- rbind(
    po[, .(dataset, x = exp(truth), what = "True CCW HR at theta*")],
    po[, .(dataset, x = exp(mean_est), what = "Average of simulated estimates (= observed)")]
  )
  d[, panel := labels[dataset]]
  lines[, panel := labels[dataset]]
  ggplot2::ggplot(d, ggplot2::aes(exp(estimate))) +
    ggplot2::geom_histogram(bins = 60, fill = "#c9dbf3", colour = "white", linewidth = 0.2) +
    ggplot2::geom_vline(
      data = lines, ggplot2::aes(xintercept = x, colour = what, linetype = what), linewidth = 0.8
    ) +
    ggplot2::facet_wrap(~panel, scales = "free", ncol = 1) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_colour_manual(values = unname(viz_cols[c("ink", "orange")]), name = NULL) +
    ggplot2::scale_linetype_manual(values = c("solid", "22"), name = NULL) +
    ggplot2::labs(x = "Estimated hazard ratio (log scale)", y = "Number of simulated studies") +
    theme_qba()
}

#' Side-by-side performance of the two datasets
#' @param perf_all Stacked performance; @param labels Dataset key -> display label.
#' @return ggplot.
plot_comparison <- function(perf_all, labels) {
  po <- perf_all[scenario == "observed" & dataset %in% names(labels)]
  pn <- perf_all[scenario == "null" & dataset %in% names(labels)]
  d <- rbind(
    po[, .(dataset, metric = "Bias (% of SE), calibrated effect",
           value = 100 * std_bias, lo = 100 * (bias - 1.96 * bias_mcse_total) / emp_se,
           hi = 100 * (bias + 1.96 * bias_mcse_total) / emp_se, ref = 0)],
    po[, .(dataset, metric = "95% CI coverage (%), calibrated effect",
           value = 100 * coverage, lo = 100 * (coverage - 1.96 * coverage_mcse),
           hi = 100 * (coverage + 1.96 * coverage_mcse), ref = 95)],
    pn[, .(dataset, metric = "Type I error (%), no effect",
           value = 100 * reject_0, lo = 100 * (reject_0 - 1.96 * reject_0_mcse),
           hi = 100 * (reject_0 + 1.96 * reject_0_mcse), ref = 5)]
  )
  d[, label := factor(labels[dataset], levels = rev(labels))]
  d[, metric := factor(metric, levels = unique(metric))]
  ggplot2::ggplot(d, ggplot2::aes(value, label)) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = ref), linetype = "22", colour = viz_cols[["ink2"]]) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = hi), height = 0, colour = viz_cols[["blue"]],
                            linewidth = 0.7) +
    ggplot2::geom_point(colour = viz_cols[["blue"]], size = 2.8) +
    ggplot2::geom_text(ggplot2::aes(label = formatC(value, format = "f", digits = 1)),
                       vjust = -1.1, size = 3.2, colour = viz_cols[["ink"]]) +
    ggplot2::facet_wrap(~metric, scales = "free_x", ncol = 1) +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_qba()
}

#' Save a ggplot and return the file path (format = "file" target)
#' @param plot ggplot; @param file Output path; @param width,height Inches.
#' @return File path.
save_plot <- function(plot, file, width = 7, height = 4.2) {
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(file, plot, width = width, height = height, dpi = 150, bg = "white")
  file
}

#' Example patients for the protocol diagram (common data format)
#' @param name Dataset key.
#' @return data.table: label, s, time, status, trt_time.
protocol_examples <- function(name) {
  d <- switch(name,
    maringe = data.table::data.table(
      label = c(
        "A  Surgery day 30, alive at 1 year", "B  Surgery day 1, dies day 250",
        "C  Surgery day 120, dies day 300", "D  No surgery, dies day 100",
        "E  No surgery, alive at 1 year"
      ),
      s = c(30, 1, 120, Inf, Inf), time = c(365.24, 250, 300, 100, 365.24), status = c(0, 1, 1, 1, 0)
    ),
    synpuf = data.table::data.table(
      label = c(
        "A  Starts day 10, event day 120", "B  Starts day 27, no event",
        "C  Starts day 60, no event", "D  Starts day 86, event day 150",
        "E  Never starts, no event", "F  Never starts, event day 15"
      ),
      s = c(10, 27, 60, 86, Inf, Inf), time = c(120, 180, 180, 150, 180, 15), status = c(1, 0, 0, 1, 0, 1)
    )
  )
  d[, `:=`(id = label, trt_time = data.table::fifelse(s <= time, s, Inf))]
  d[]
}

#' Clone segments for the protocol diagram, using the pipeline's censoring rules
#' @param ex Example patients; @param name Dataset key; @param spec Pipeline spec.
#' @return list(seg = segments with weight status, ends = end markers, arms = arm labels).
protocol_segments <- function(ex, name, spec) {
  if (name == "maringe") {
    g <- spec$grace
    cl <- build_clones(ex, g)
    cl[, arm_lab := data.table::fifelse(arm == 1L, "Surgery within 6 months", "No surgery within 6 months")]
    seg <- rbind(
      # Control arm: weight 1 / P(not yet operated | X), updated at every surgery date
      cl[arm == 0L, .(label, arm_lab, x0 = 0, x1 = end, wt = "Up-weighted (IPCW)")],
      cl[arm == 1L, .(label, arm_lab, x0 = 0, x1 = pmin(end, g), wt = "Weight 1")],
      cl[arm == 1L & end > g, .(label, arm_lab, x0 = g, x1 = end, wt = "Up-weighted (IPCW)")]
    )
    win <- data.table::data.table(arm_lab = unique(cl$arm_lab), w = g)
  } else {
    wins <- unlist(spec$windows)
    cl <- data.table::rbindlist(lapply(wins, function(w) {
      x <- data.table::copy(ex)
      # Same rule as ccw_window_arm(): censored at w if followed beyond w without starting
      x[, nostart := as.integer(time > w & trt_time > w)]
      x[, `:=`(
        end = data.table::fifelse(nostart == 1L, w, time),
        event = data.table::fifelse(nostart == 1L, 0L, as.integer(status)),
        art_cens = nostart, w = w,
        arm_lab = paste0("Start within ", w, " days", if (w == wins[1]) " (reference)" else "")
      )]
    }))
    seg <- rbind(
      cl[, .(label, arm_lab, x0 = 0, x1 = pmin(end, w), wt = "Weight 1")],
      cl[end > w, .(label, arm_lab, x0 = w, x1 = end,
                    wt = data.table::fifelse(trt_time >= w - spec$recent_days,
                                             "Up-weighted (IPCW)", "Weight 1"))]
    )
    win <- unique(cl[, .(arm_lab, w)])
  }
  ends <- cl[, .(label, arm_lab, x = end, what = data.table::fcase(
    event == 1L, "Outcome event",
    art_cens == 1L, "Artificially censored",
    default = "End of follow-up"
  ))]
  starts <- unique(cl[is.finite(trt_time), .(label, arm_lab, x = trt_time)])
  # Strategy of interest first, reference strategy second
  arms <- rev(sort(unique(cl$arm_lab)))
  if (name == "synpuf") arms <- rev(arms)
  list(seg = seg, ends = ends, starts = starts, win = win, arms = arms)
}

#' Protocol diagram: how example patients are cloned, censored and weighted
#' @param name Dataset key; @param spec Pipeline spec; @param horizon End of follow-up (days).
#' @return ggplot.
plot_protocol <- function(name, spec, horizon) {
  ex <- protocol_examples(name)
  ps <- protocol_segments(ex, name, spec)
  lev <- rev(ex$label)
  for (nm in c("seg", "ends", "starts")) ps[[nm]][, label := factor(label, levels = lev)]
  for (nm in c("seg", "ends", "starts", "win")) ps[[nm]][, arm_lab := factor(arm_lab, levels = ps$arms)]
  ggplot2::ggplot() +
    ggplot2::geom_rect(data = ps$win, ggplot2::aes(xmin = 0, xmax = w, ymin = -Inf, ymax = Inf),
                       fill = "#efeee9") +
    ggplot2::geom_segment(data = ps$seg, ggplot2::aes(x = x0, xend = x1, y = label, yend = label, colour = wt),
                          linewidth = 2.6) +
    ggplot2::geom_point(data = ps$starts, ggplot2::aes(x = x, y = label, shape = "Treatment start"),
                        size = 2.6, colour = viz_cols[["ink"]], fill = "white") +
    ggplot2::geom_point(data = ps$ends, ggplot2::aes(x = x, y = label, shape = what, size = what),
                        colour = viz_cols[["ink"]], stroke = 1.1) +
    ggplot2::scale_size_manual(
      values = c("Outcome event" = 3, "Artificially censored" = 7, "End of follow-up" = 3),
      guide = "none"
    ) +
    ggplot2::facet_wrap(~arm_lab, ncol = 2) +
    ggplot2::scale_colour_manual(
      values = c("Weight 1" = viz_cols[["blue"]], "Up-weighted (IPCW)" = viz_cols[["orange"]]),
      name = NULL
    ) +
    ggplot2::scale_shape_manual(
      values = c("Treatment start" = 24, "Outcome event" = 4, "Artificially censored" = 124,
                 "End of follow-up" = 15),
      name = NULL
    ) +
    ggplot2::scale_x_continuous(limits = c(0, horizon * 1.02), expand = c(0, 0)) +
    ggplot2::labs(x = "Days since time zero", y = NULL) +
    # Two rows so the symbol legend fits the figure width
    ggplot2::guides(shape = ggplot2::guide_legend(nrow = 2)) +
    theme_qba() +
    ggplot2::theme(
      text = ggplot2::element_text(size = 13), axis.text = ggplot2::element_text(size = 11),
      panel.grid.major.y = ggplot2::element_blank(), legend.box = "vertical"
    )
}
