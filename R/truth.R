# True value psi(theta): the large-sample limit of the CCW estimator in the
# simulated world (the "CCW model parameter"; ADEMP E1), and the calibration of theta*.
# One large dataset per (theta, seed) is the unit of work, so truths and calibration steps
# run as parallel branches, and every evaluation is reused wherever the same (theta, seed)
# is needed again (CRN: seed j is the same large dataset at every theta).

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
#' @return One-row data.table: theta, seed_id, estimate, n_truth.
truth_point <- function(dgm, data, theta, seed_id, seeds, pipeline, cfg, known = NULL) {
  checkmate::assert_number(theta)
  checkmate::assert_int(seed_id, lower = 1, upper = length(seeds))
  if (!is.null(known) && nrow(known)) {
    # Local copies: inside data.table's [ ] the bare names would refer to the columns
    th <- theta
    sid <- seed_id
    hit <- known[abs(theta - th) < 1e-10 & seed_id == sid]
    if (nrow(hit)) return(hit[1, .(theta, seed_id, estimate, n_truth)])
  }
  big <- stack_copies(data, ceiling(cfg$truth$n_target / nrow(data)))
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
#' @param dgm,data,theta,pipeline,cfg As in `truth_point()`; @param n_rep Large datasets.
#' @return One-row data.table: theta, truth, truth_mcse, n_truth, n_rep.
compute_truth <- function(dgm, data, theta, pipeline, cfg, n_rep = cfg$truth$n_rep) {
  # Same seeds at every theta (CRN), so psi(theta) is smooth for calibration
  seeds <- make_crn_seeds(cfg$seed_truth, n_rep)
  rows <- data.table::rbindlist(lapply(seq_len(n_rep), function(j) {
    truth_point(dgm, data, theta, j, seeds, pipeline, cfg)
  }))
  summarise_truth(rows)
}

#' Starting values of the calibration search (paper Step 3)
#' @param dgm Fitted DGM; @param min_gap Smallest |theta_hat| for which 0 is the second start.
#' @return Two theta values.
calib_start <- function(dgm, min_gap = 0.25) {
  th <- dgm$theta_hat
  # psi(0) is needed anyway as the null-scenario truth, so 0 is a free second start; it is
  # only uninformative for the secant when theta_hat is itself close to 0
  if (abs(th) >= min_gap) c(th, 0) else c(th, th - 0.5)
}

#' Next calibration point from all evaluations so far
#' @param rows From `truth_point()`; @param target Observed CCW estimate.
#' @param cc `cfg$calibration` (tol, bound).
#' @return Scalar theta: the best point so far once within tol, else a secant step.
calib_next <- function(rows, target, cc) {
  checkmate::assert_number(target, finite = TRUE)
  checkmate::assert_number(cc$tol, lower = 0)
  checkmate::assert_number(cc$bound, lower = 0)
  ps <- summarise_truth(rows)
  ps[, gap := truth - target]
  best <- ps[which.min(abs(gap))]
  if (abs(best$gap) < cc$tol || nrow(ps) < 2) return(best$theta)
  # Interpolate between the closest points bracketing the target, else secant through the
  # two closest points; psi is smooth and monotone in theta under CRN
  lo <- ps[gap < 0][which.max(gap)]
  hi <- ps[gap > 0][which.min(gap)]
  pair <- if (nrow(lo) && nrow(hi)) rbind(lo, hi) else ps[order(abs(gap))][1:2]
  slope <- diff(pair$truth) / diff(pair$theta)
  th <- pair$theta[1] + (target - pair$truth[1]) / slope
  if (!is.finite(th)) return(best$theta)
  max(-cc$bound, min(cc$bound, th))
}

#' Finish the calibration: serial secant steps if the parallel stages did not converge
#' @param rows Evaluations from the parallel stages; @param dgm,data,pipeline,cfg As above.
#' @param obs_result Observed pipeline result; @param seeds Truth seeds.
#' @param n_stages Parallel refinement stages already run (count towards `max_iter`).
#' @return List: theta_star, psi_star, at_bound, path (search history), evals (all rows).
finish_calibration <- function(rows, dgm, data, obs_result, pipeline, cfg, seeds, n_stages = 2L) {
  cc <- cfg$calibration
  target <- obs_result$estimate
  ids <- seq_len(cc$n_rep)
  for (i in seq_len(max(0L, cc$max_iter - n_stages))) {
    th <- calib_next(rows, target, cc)
    if (any(abs(rows$theta - th) < 1e-10)) break
    rows <- rbind(rows, data.table::rbindlist(lapply(ids, function(j) {
      truth_point(dgm, data, th, j, seeds, pipeline, cfg)
    })))
  }
  ps <- summarise_truth(rows, ids)
  best <- ps[which.min(abs(truth - target))]
  # Order of first evaluation = search order
  first <- unique(rows[, .(theta)])
  path <- ps[match(first$theta, ps$theta), .(step = seq_len(.N), theta, psi = truth)]
  list(
    theta_star = best$theta,
    psi_star = best$truth,
    at_bound = abs(best$theta) >= cc$bound,
    path = path,
    evals = rows
  )
}

#' Map target psi values to theta through the evaluated psi(theta) points
#' @param psi Target values; @param pts data.table theta, psi (at least two distinct thetas).
#' @param bound |theta| cap.
#' @return theta values (linear interpolation; linear extrapolation with the fitted slope).
psi_to_theta <- function(psi, pts, bound = Inf) {
  pts <- unique(pts[, .(theta, psi)], by = "theta")[order(psi)]
  checkmate::assert_true(nrow(pts) >= 2)
  slope <- unname(stats::coef(stats::lm(psi ~ theta, data = pts))[2])
  th <- stats::approx(pts$psi, pts$theta, xout = psi, rule = 2)$y
  below <- psi < pts$psi[1]
  above <- psi > pts$psi[nrow(pts)]
  th[below] <- pts$theta[1] + (psi[below] - pts$psi[1]) / slope
  th[above] <- pts$theta[nrow(pts)] + (psi[above] - pts$psi[nrow(pts)]) / slope
  pmax(-bound, pmin(bound, th))
}
