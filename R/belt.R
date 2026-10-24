# Confidence belt: CIs by inverting a test over a grid of true effects (Neyman construction).
# Unlike the bootstrap-t, the quantiles of T may change with the true effect. The grid is placed
# where CI endpoints fall (Chebyshev-Lobatto nodes in psi), and replicates per grid point are
# added only until the quantiles reach a target MCSE.

#' Range of psi the belt must cover: where CI endpoints of the observed study and of the
#' checked replicates fall
#' @param obs_result Observed pipeline result; @param sim_results Main replicates (Wald CIs).
#' @param tail Share of replicate endpoints allowed beyond the range, per side.
#' @param expand Widening of the range on each side, as a share of its width.
#' @return c(lo, hi).
belt_psi_range <- function(obs_result, sim_results, tail, expand) {
  d <- sim_results[is.finite(ci_lo) & is.finite(ci_hi)]
  # Wald endpoints locate the belt endpoints closely; `expand` absorbs the difference
  lo <- min(obs_result$ci_lo, stats::quantile(d$ci_lo, tail, names = FALSE))
  hi <- max(obs_result$ci_hi, stats::quantile(d$ci_hi, 1 - tail, names = FALSE))
  c(lo, hi) + c(-1, 1) * expand * (hi - lo)
}

#' Chebyshev-Lobatto nodes on an interval (endpoints included)
#' @param range c(lo, hi); @param k Number of nodes (>= 2).
#' @return Increasing numeric vector.
cheb_lobatto <- function(range, k) {
  checkmate::assert_int(k, lower = 2)
  # Denser near the ends, where a polynomial fit is least stable
  mean(range) - diff(range) / 2 * cos(pi * (0:(k - 1)) / (k - 1))
}

#' Theta values of the extra belt grid points
#' @param range From `belt_psi_range()`; @param truth_table Main-scenario truths.
#' @param psi_pts Large-dataset rows from `truth_point()` (calibration and main truths).
#' @param cfg_belt Belt config (n_nodes, min_gap, bound).
#' @return Sorted theta values, excluding nodes the main scenarios already cover.
belt_theta_plan <- function(range, truth_table, psi_pts, cfg_belt) {
  nodes <- cheb_lobatto(range, cfg_belt$n_nodes)
  # The main scenarios (10,000 replicates each) stand in for nodes close to them
  keep <- vapply(nodes, function(x) all(abs(x - truth_table$truth) > cfg_belt$min_gap),
                 logical(1))
  pts <- summarise_truth(psi_pts)[, .(theta, psi = truth)]
  th <- unique(round(psi_to_theta(nodes[keep], pts, cfg_belt$bound), 6))
  sort(th[vapply(th, function(x) all(abs(x - truth_table$theta) > 1e-6), logical(1))])
}

#' Monte Carlo SE of a sample quantile (distribution-free, from order statistics)
#' @param x Sample; @param p Probability.
#' @return Scalar.
quantile_mcse <- function(x, p) {
  x <- sort(x[is.finite(x)])
  n <- length(x)
  z <- stats::qnorm(0.975)
  h <- z * sqrt(n * p * (1 - p))
  j <- pmin(n, pmax(1, round(n * p + c(-h, h))))
  (x[j[2]] - x[j[1]]) / (2 * z)
}

#' Quantiles of T = (estimate - truth) / se at each grid point
#' @param d Replicates with columns theta, rep, estimate, se, truth.
#' @param level Nominal level.
#' @return data.table: theta, psi, n, q_lo, q_hi, q_lo_mcse, q_hi_mcse.
belt_quantiles <- function(d, level) {
  a <- (1 - level) / 2
  d[, .(
    n = .N,
    q_lo = stats::quantile((estimate - truth) / se, a, names = FALSE, type = 8),
    q_hi = stats::quantile((estimate - truth) / se, 1 - a, names = FALSE, type = 8),
    q_lo_mcse = quantile_mcse((estimate - truth) / se, a),
    q_hi_mcse = quantile_mcse((estimate - truth) / se, 1 - a)
  ), by = .(theta, psi = truth)][order(psi)]
}

#' Replicates needed per grid point for the target quantile MCSE (adaptive B)
#' @param sims First-stage belt replicates; @param truths data.table theta, truth.
#' @param cfg_belt Belt config (B_start, B_max, target_mcse).
#' @param batch_size Replicates per batch.
#' @param level Nominal level.
#' @return data.table: theta, n, mcse, n_batches (total batches to run at that theta).
belt_topup_plan <- function(sims, truths, cfg_belt, batch_size, level = 0.95) {
  d <- merge(sims[!is.na(estimate) & se > 0, .(theta, rep, estimate, se)],
             truths[, .(theta, truth)], by = "theta")
  q <- belt_quantiles(d, level)
  q[, mcse := pmax(q_lo_mcse, q_hi_mcse)]
  # MCSE of a quantile falls as 1 / sqrt(n)
  q[, n_req := ceiling(n * (mcse / cfg_belt$target_mcse)^2)]
  q[, n_req := pmin(cfg_belt$B_max, pmax(cfg_belt$B_start, n_req))]
  q[, .(theta, n, mcse, n_batches = as.integer(ceiling(n_req / batch_size)))]
}

#' Run one top-up batch at one grid point, or nothing if the plan does not need it
#' @param data,dgm,theta,seeds,pipeline,batch_size As in `run_batch()`.
#' @param batch_id Batch index (continues after the first stage).
#' @param plan From `belt_topup_plan()`.
#' @return Replicate rows (possibly none).
run_belt_topup <- function(data, dgm, theta, seeds, batch_id, plan, pipeline, batch_size) {
  need <- plan$n_batches[abs(plan$theta - theta) < 1e-10]
  if (!length(need) || batch_id > need) {
    return(data.table::data.table(theta = numeric(), rep = integer(), estimate = numeric(),
                                  se = numeric()))
  }
  run_batch(data, dgm, theta, seeds, batch_id, pipeline, batch_size)
}

#' Smooth quantile curves q(psi): weighted polynomial, constant beyond the grid
#' @param q From `belt_quantiles()`; @param degree Polynomial degree.
#' @return list(lo, hi = coefficients, range, degree).
belt_curve <- function(q, degree) {
  deg <- min(degree, nrow(q) - 1L)
  X <- outer(q$psi, 0:deg, `^`)
  fit <- function(y) stats::lm.wfit(X, y, w = q$n)$coefficients
  list(lo = fit(q$q_lo), hi = fit(q$q_hi), range = range(q$psi), degree = deg)
}

#' Evaluate a quantile curve
#' @param curve From `belt_curve()`; @param which "lo" or "hi"; @param psi Values.
#' @return Numeric vector.
belt_q <- function(curve, which, psi) {
  # No extrapolation of the polynomial: outside the grid the nearest quantile is used
  p <- pmin(pmax(psi, curve$range[1]), curve$range[2])
  drop(outer(p, 0:curve$degree, `^`) %*% curve[[which]])
}

#' Invert the belt: all psi0 with q_lo(psi0) <= (estimate - psi0) / se <= q_hi(psi0)
#' @param estimate,se Vectors; @param curve From `belt_curve()`.
#' @return list(lo, hi, outside = endpoint beyond the simulated psi range).
belt_invert <- function(estimate, se, curve) {
  # Vectorised bisection; (estimate - psi) / se - q(psi) decreases in psi when the curves
  # are flatter than 1 / se, which holds by a wide margin here
  root <- function(which) {
    a <- estimate - 10 * se
    b <- estimate + 10 * se
    for (i in 1:60) {
      m <- (a + b) / 2
      pos <- (estimate - m) / se - belt_q(curve, which, m) > 0
      a <- ifelse(pos, m, a)
      b <- ifelse(pos, b, m)
    }
    (a + b) / 2
  }
  lo <- root("hi")
  hi <- root("lo")
  list(lo = lo, hi = hi, outside = lo < curve$range[1] | hi > curve$range[2])
}

#' Build the belt: grid quantiles and curves, on all replicates and per cross-fitting fold
#' @param sims Replicates at theta*, 0 and the extra grid; @param truths data.table theta, truth.
#' @param cfg_belt Belt config; @param level Nominal level.
#' @return list(quantiles, curve (all replicates), data (replicates, for the coverage check),
#'   degree, level).
fit_belt <- function(sims, truths, cfg_belt, level = 0.95) {
  d <- merge(sims[!is.na(estimate) & se > 0, .(theta, rep, estimate, se)],
             truths[, .(theta, truth)], by = "theta")
  q_all <- belt_quantiles(d, level)
  list(quantiles = q_all, curve = belt_curve(q_all, cfg_belt$degree), data = d,
       degree = cfg_belt$degree, level = level)
}

#' Coverage of the belt CI in the main scenarios, with the checked grid point left out
#' @param sim_results Main replicates; @param truth_table From `make_truth_table()`.
#' @param belt From `fit_belt()`.
#' @return Rows in the `compare_ci_methods()` format, plus share_outside.
belt_coverage <- function(sim_results, truth_table, belt) {
  d <- merge(sim_results[!is.na(estimate) & se > 0], truth_table[, .(dataset, theta, scenario, truth)],
             by = c("dataset", "theta"))
  data.table::rbindlist(lapply(c("observed", "null"), function(tg) {
    th <- d[scenario == tg, theta[1]]
    res <- data.table::rbindlist(lapply(0:1, function(k) {
      # A CI covers psi0 exactly when the test at psi0 does not reject, so with psi0 on the grid
      # coverage holds by construction. Dropping that grid point makes the belt interpolate, the
      # situation of a real study. Folds split replicates (CRN links replicate b across thetas).
      g <- belt$data[abs(theta - th) > 1e-8 & rep %% 2L != k]
      cv <- belt_curve(belt_quantiles(g, belt$level), belt$degree)
      x <- d[scenario == tg & rep %% 2L == k]
      ci <- belt_invert(x$estimate, x$se, cv)
      data.table::data.table(lo = ci$lo, hi = ci$hi, truth = x$truth, outside = ci$outside)
    }))
    hit <- res$lo <= res$truth & res$truth <= res$hi
    m <- mean(hit)
    data.table::data.table(
      method = "Confidence belt", scenario = tg, source = "grid", check = "leave-one-out",
      coverage = m, mcse = sqrt(m * (1 - m) / nrow(res)),
      miss_below = mean(res$hi < res$truth), miss_above = mean(res$lo > res$truth),
      n = nrow(res), share_outside = mean(res$outside)
    )
  }))
}

#' Figure: quantiles of T across the grid with the fitted curves
#' @param belt From `fit_belt()`; @param level Nominal level.
#' @return ggplot.
plot_belt <- function(belt, level = 0.95) {
  q <- data.table::melt(belt$quantiles, id.vars = c("theta", "psi", "n"),
                        measure.vars = c("q_lo", "q_hi"), variable.name = "which", value.name = "q")
  grid <- seq(belt$curve$range[1], belt$curve$range[2], length.out = 200)
  cv <- data.table::rbindlist(lapply(c("lo", "hi"), function(w) {
    data.table::data.table(which = paste0("q_", w), psi = grid, q = belt_q(belt$curve, w, grid))
  }))
  zz <- stats::qnorm(1 - (1 - level) / 2)
  lab <- c(q_lo = "2.5% quantile", q_hi = "97.5% quantile")
  ggplot2::ggplot(q, ggplot2::aes(exp(psi), q, colour = which)) +
    ggplot2::geom_hline(yintercept = c(-zz, zz), linetype = "dashed") +
    ggplot2::geom_line(data = cv, linewidth = 0.9) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_colour_manual(values = c(q_lo = "#1f6fb4", q_hi = "red"), labels = lab, name = NULL) +
    ggplot2::labs(x = "True CCW HR (grid point)", y = "Quantile of T = (estimate - truth) / SE") +
    theme_paper(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}

#' Save the belt figure (PNG and EPS)
#' @param belt From `fit_belt()`; @param dir Output folder.
#' @return File paths.
save_belt_figure <- function(belt, dir) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  p <- plot_belt(belt)
  png <- file.path(dir, "fig_belt.png")
  eps <- file.path(dir, "fig_belt.eps")
  ggplot2::ggsave(png, p, width = 7, height = 5, dpi = 600, bg = "white")
  ggplot2::ggsave(eps, p, width = 7, height = 5, device = grDevices::cairo_ps)
  c(png, eps)
}
