# True value psi(theta): the large-sample limit of the CCW estimator in the simulated world
# (ADEMP E1), and the inverse calibration of the generating value (paper Section 2.7.1).

#' Stack K copies of the base data with unique ids
#' @param data Base data; @param K Number of copies.
#' @return data.table with K * nrow(data) rows.
stack_copies <- function(data, K) {
  big <- data[rep(seq_len(.N), K)]
  big[, id := paste0(id, "_", rep(seq_len(K), each = nrow(data)))]
  big[]
}

#' Large-sample truth at one theta
#' @param dgm Fitted DGM; @param data Base data; @param theta Treatment log-HR in the DGM.
#' @param pipeline Pipeline spec; @param cfg Full config.
#' @param n_rep,n_target Large datasets to average and patients per dataset.
#' @param seed Seed of the large datasets (CRN across theta).
#' @return One-row data.table: theta, truth, truth_mcse, n_truth, n_rep.
compute_truth <- function(dgm, data, theta, pipeline, cfg, n_rep = cfg$truth$n_rep,
                          n_target = cfg$truth$n_target, seed = cfg$seed_truth) {
  checkmate::assert_number(theta)
  big <- stack_copies(data, ceiling(n_target / nrow(data)))
  # Same seeds at every theta (CRN), so psi(theta) is smooth in theta
  seeds <- make_crn_seeds(seed, n_rep)
  est <- vapply(seeds, function(sd) {
    rng <- draw_replicate_rng(sd, nrow(big))
    sim <- simulate_outcomes(big, dgm, theta, rng$u_event, rng$u_cens)
    run_pipeline(sim, pipeline, robust = FALSE)$estimate
  }, numeric(1))
  data.table::data.table(
    theta = theta,
    truth = mean(est),
    truth_mcse = if (n_rep > 1) stats::sd(est) / sqrt(n_rep) else NA_real_,
    n_truth = nrow(big),
    n_rep = n_rep
  )
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
