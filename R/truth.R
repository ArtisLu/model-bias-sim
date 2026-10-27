# True value psi(theta): the large-sample limit of the CCW estimator in the simulated world
# (ADEMP E1), and the inverse calibration of the generating value (paper Section 2.7.1).
# One large dataset per (theta, seed) is the unit of work, so truths run as parallel branches and
# every evaluation is reused wherever the same (theta, seed) is needed again (CRN: seed j is the
# same large dataset at every theta).

#' Stack K copies of the base data with unique ids
#' @param data Base data; @param K Number of copies.
#' @return data.table with K * nrow(data) rows.
stack_copies <- function(data, K) {
  big <- data[rep(seq_len(.N), K)]
  big[, id := paste0(id, "_", rep(seq_len(K), each = nrow(data)))]
  big[]
}

#' CCW estimate on one large dataset (one truth seed at one theta)
#' @param dgm Fitted DGM; @param data Base data; @param theta Treatment log-HR in the DGM.
#' @param seed_id Index into `seeds`; @param seeds Truth seeds from `make_crn_seeds()`.
#' @param pipeline Pipeline spec; @param cfg Full config (uses `truth$n_target`).
#' @param known Earlier evaluations (theta, seed_id, estimate, n_truth) to reuse, or NULL.
#' @param n_target Patients per large dataset (default `cfg$truth$n_target`).
#' @return One-row data.table: theta, seed_id, estimate, n_truth.
truth_point <- function(dgm, data, theta, seed_id, seeds, pipeline, cfg, known = NULL,
                        n_target = cfg$truth$n_target) {
  checkmate::assert_number(theta)
  checkmate::assert_int(seed_id, lower = 1, upper = length(seeds))
  if (!is.null(known) && nrow(known)) {
    # Local copies: inside data.table's [ ] the bare names would refer to the columns
    th <- theta
    sid <- seed_id
    hit <- known[abs(theta - th) < 1e-10 & seed_id == sid]
    if (nrow(hit)) return(hit[1, .(theta, seed_id, estimate, n_truth)])
  }
  big <- stack_copies(data, ceiling(n_target / nrow(data)))
  rng <- draw_replicate_rng(seeds[seed_id], nrow(big))
  sim <- simulate_outcomes(big, dgm, theta, rng$u_event, rng$u_cens)
  data.table::data.table(
    theta = theta, seed_id = seed_id,
    estimate = run_pipeline(sim, pipeline, robust = FALSE)$estimate, n_truth = nrow(big)
  )
}

#' Average large-dataset estimates into true values
#' @param rows From `truth_point()`; @param seed_ids Seeds to use (default: all present).
#' @return data.table, one row per theta: theta, truth, truth_mcse, n_truth, n_rep.
summarise_truth <- function(rows, seed_ids = NULL) {
  if (!is.null(seed_ids)) rows <- rows[seed_id %in% seed_ids]
  rows <- unique(rows, by = c("theta", "seed_id"))
  rows[, .(
    truth = mean(estimate),
    truth_mcse = if (.N > 1) stats::sd(estimate) / sqrt(.N) else NA_real_,
    n_truth = n_truth[1], n_rep = .N
  ), by = theta][order(theta)]
}

#' Large-sample truth at one theta (serial; the pipeline branches over seeds instead)
#' @param dgm,data,theta,pipeline,cfg As in `truth_point()`; @param n_rep Large datasets to average.
#' @param n_target Patients per large dataset; @param seed Seed of the large datasets (CRN across theta).
#' @return One-row data.table: theta, truth, truth_mcse, n_truth, n_rep.
compute_truth <- function(dgm, data, theta, pipeline, cfg, n_rep = cfg$truth$n_rep,
                          n_target = cfg$truth$n_target, seed = cfg$seed_truth) {
  # Same seeds at every theta (CRN), so psi(theta) is smooth in theta
  seeds <- make_crn_seeds(seed, n_rep)
  rows <- data.table::rbindlist(lapply(seq_len(n_rep), function(j) {
    truth_point(dgm, data, theta, j, seeds, pipeline, cfg, n_target = n_target)
  }))
  summarise_truth(rows)
}

#' Pipeline estimates of simulated studies at one theta (failed fits are NA)
#' @param data Base data; @param dgm Fitted DGM; @param theta Generating log-HR.
#' @param seeds Replicate seeds; @param pipeline Pipeline spec.
#' @return Numeric vector, one estimate per seed.
replicate_estimates <- function(data, dgm, theta, seeds, pipeline) {
  vapply(seeds, function(sd) {
    rng <- draw_replicate_rng(sd, nrow(data))
    sim <- simulate_outcomes(data, dgm, theta, rng$u_event, rng$u_cens)
    # Only the point estimate is needed, which does not depend on the variance estimator
    safely_run(run_pipeline, sim, pipeline, robust = FALSE)$estimate
  }, numeric(1))
}

#' Inverse calibration: theta whose mean simulated estimate m(theta) equals a target (paper 2.7.1)
#' The search starts at `start[1]` and stops there if m is already within `cc$tol` of the target.
#' @param dgm Fitted DGM; @param data Base data; @param target Estimate to reproduce.
#' @param seeds Replicate seeds, reused at every theta (CRN makes m(theta) smooth and monotone).
#' @param pipeline Pipeline spec; @param cc List with max_iter, tol, bound.
#' @param start First theta (default `cc$start`); an optional second element is the second point.
#' @param slope Expected dm/dtheta; if given, the second point is a Newton step from `start[1]`;
#'   without it the second point is `start[2]` or, if absent, a unit-slope step.
#' @return List: theta_cal, mean_est, mean_est_mcse, n_ok, estimates (at theta_cal), slope,
#'   at_bound, n_eval, path.
calibrate_mean <- function(dgm, data, target, seeds, pipeline, cc, start = cc$start, slope = NULL) {
  checkmate::assert_number(target)
  checkmate::assert_numeric(start, min.len = 1, any.missing = FALSE)
  clamp <- function(x) max(-cc$bound, min(cc$bound, x))
  evals <- list()
  run_at <- function(th) {
    est <- replicate_estimates(data, dgm, th, seeds, pipeline)
    evals[[length(evals) + 1L]] <<- list(theta = th, m = mean(est, na.rm = TRUE), est = est)
  }
  prior_slope <- if (is.null(slope)) 1 else slope
  last_slope <- function() {
    n <- length(evals)
    if (n < 2L) return(prior_slope)
    s <- (evals[[n]]$m - evals[[n - 1L]]$m) / (evals[[n]]$theta - evals[[n - 1L]]$theta)
    # A flat or undefined secant (failed fits, same theta) falls back to the prior slope
    if (is.finite(s) && abs(s) > 1e-8) s else prior_slope
  }
  run_at(clamp(start[1]))
  if (is.finite(evals[[1]]$m) && abs(evals[[1]]$m - target) >= cc$tol) {
    th2 <- if (!is.null(slope)) {
      clamp(evals[[1]]$theta + (target - evals[[1]]$m) / slope)
    } else if (length(start) >= 2L) {
      start[2]
    } else {
      clamp(evals[[1]]$theta + (target - evals[[1]]$m))
    }
    if (abs(th2 - evals[[1]]$theta) > 1e-6) run_at(th2)
  }
  for (i in seq_len(cc$max_iter)) {
    n <- length(evals)
    if (!is.finite(evals[[n]]$m) || abs(evals[[n]]$m - target) < cc$tol) break
    th_new <- clamp(evals[[n]]$theta + (target - evals[[n]]$m) / last_slope())
    if (abs(th_new - evals[[n]]$theta) < 1e-6) break
    run_at(th_new)
  }
  path <- data.table::rbindlist(lapply(seq_along(evals), function(k) {
    data.table::data.table(step = k, theta = evals[[k]]$theta, mean_est = evals[[k]]$m,
                           n_ok = sum(!is.na(evals[[k]]$est)))
  }))
  if (!any(is.finite(path$mean_est))) {
    cli::cli_abort("All simulated analyses failed in calibration.")
  }
  best <- which.min(abs(path$mean_est - target))
  est <- evals[[best]]$est
  list(
    theta_cal = path$theta[best],
    mean_est = path$mean_est[best],
    mean_est_mcse = stats::sd(est, na.rm = TRUE) / sqrt(sum(!is.na(est))),
    n_ok = sum(!is.na(est)),
    estimates = est,
    slope = last_slope(),
    at_bound = abs(path$theta[best]) >= cc$bound,
    n_eval = nrow(path),
    path = path
  )
}

#' Calibrate the observed analysis: m(theta) = observed estimate (paper Steps 2.3 and 2.7.1)
#' @param dgm Fitted DGM; @param data Base data; @param obs_result Observed pipeline result.
#' @param seeds CRN seeds of the main run (the first `cfg$calibration$B` are used).
#' @param pipeline Pipeline spec; @param cfg Full config.
#' @param slope,tol_scale Dataset-level guess of dm/dtheta (first Newton step) and the multiplier
#'   of the tolerance for the scale of the estimate (`calib_slope`, `tol_scale` in the config).
#' @return `calibrate_mean()` result without the replicate estimates, plus target.
calibrate_observed <- function(dgm, data, obs_result, seeds, pipeline, cfg, slope = NULL,
                               tol_scale = 1) {
  cc <- cfg$calibration
  cc$tol <- cc$tol * tol_scale
  cal <- calibrate_mean(dgm, data, obs_result$estimate, seeds[seq_len(min(cc$B, length(seeds)))],
                        pipeline, cc, slope = slope)
  cal$estimates <- NULL
  c(cal, list(target = obs_result$estimate, B = min(cc$B, length(seeds))))
}

#' Monotone interpolation of psi(theta) through the evaluated points
#' @param pts data.table theta, psi (at least two distinct thetas).
#' @return Function of theta (Hyman monotone cubic spline; flat beyond the points).
psi_curve <- function(pts) {
  pts <- unique(pts[, .(theta, psi)], by = "theta")[order(theta)]
  checkmate::assert_true(nrow(pts) >= 2)
  # CRN keeps psi monotone in theta; cummax only guards against Monte Carlo ties
  y <- cummax(pts$psi)
  stats::splinefun(pts$theta, y, method = "hyman")
}

#' Map target psi values to theta by inverting `psi_curve()`
#' @param psi Target values; @param f From `psi_curve()`; @param th_range Theta range searched.
#' @return theta values; targets beyond psi at the range ends map to the range ends.
psi_to_theta <- function(psi, f, th_range) {
  ends <- f(th_range)
  vapply(psi, function(p) {
    if (p <= ends[1]) return(th_range[1])
    if (p >= ends[2]) return(th_range[2])
    stats::uniroot(function(t) f(t) - p, th_range, tol = 1e-9)$root
  }, numeric(1))
}
