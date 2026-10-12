# True value psi(theta): the large-sample limit of the CCW estimator in the
# simulated world (the "CCW model parameter"; ADEMP E1), and the calibration of theta*.

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
#' @param pipeline Pipeline spec; @param cfg Full config (uses `truth`, `seed_truth`).
#' @param n_rep Number of large datasets to average.
#' @return One-row data.table: theta, truth, truth_mcse, n_truth, n_rep.
compute_truth <- function(dgm, data, theta, pipeline, cfg, n_rep = cfg$truth$n_rep) {
  checkmate::assert_number(theta)
  big <- stack_copies(data, ceiling(cfg$truth$n_target / nrow(data)))
  # Same seeds at every theta (CRN), so psi(theta) is smooth for calibration
  seeds <- make_crn_seeds(cfg$seed_truth, n_rep)
  est <- vapply(seeds, function(sd) {
    u <- draw_replicate_rng(sd, nrow(big))$u_event
    run_pipeline(simulate_outcomes(big, dgm, theta, u), pipeline, robust = FALSE)$estimate
  }, numeric(1))
  data.table::data.table(
    theta = theta,
    truth = mean(est),
    truth_mcse = if (n_rep > 1) stats::sd(est) / sqrt(n_rep) else NA_real_,
    n_truth = nrow(big),
    n_rep = n_rep
  )
}

#' Find theta* whose CCW truth equals the observed CCW estimate (paper Step 3)
#' @param dgm Fitted DGM; @param data Base data; @param obs_result Observed pipeline result.
#' @param pipeline Pipeline spec; @param cfg Full config (uses `calibration`).
#' @return List: theta_star, psi_star (search precision), path (search history).
calibrate_theta <- function(dgm, data, obs_result, pipeline, cfg) {
  cc <- cfg$calibration
  psi <- function(th) compute_truth(dgm, data, th, pipeline, cfg, n_rep = cc$n_rep)$truth
  target <- obs_result$estimate
  # Secant search started at the fitted coefficient; psi is smooth and monotone under CRN
  th <- c(dgm$theta_hat, dgm$theta_hat - 0.5)
  ps <- vapply(th, psi, numeric(1))
  for (i in seq_len(cc$max_iter)) {
    if (abs(ps[length(ps)] - target) < cc$tol) break
    n <- length(th)
    slope <- (ps[n] - ps[n - 1]) / (th[n] - th[n - 1])
    th_new <- th[n] + (target - ps[n]) / slope
    th_new <- max(-cc$bound, min(cc$bound, th_new))
    th <- c(th, th_new)
    ps <- c(ps, psi(th_new))
  }
  best <- which.min(abs(ps - target))
  list(
    theta_star = th[best],
    psi_star = ps[best],
    at_bound = abs(th[best]) >= cc$bound,
    path = data.table::data.table(step = seq_along(th), theta = th, psi = ps)
  )
}
