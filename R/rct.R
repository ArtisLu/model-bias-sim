# Target estimand for a pipeline that is not consistent (paper Section 4.6): the log HR of a
# hypothetical randomised trial in the simulated world, analysed without weights.
#
# Construction (assumption, see protocol): every patient of the base data is a member of both arms
# (randomisation balances the covariates). Arm 1 receives surgery at the patient's own surgery date
# if it fell within the grace period, otherwise at a date drawn from the observed surgery dates
# within the grace period; arm 0 has no surgery. Outcomes are predicted by the prediction model,
# with a new prediction for every copy of a patient (stacked copies, as for the CCW truth).
# Patients without surgery up to the end of the grace period (all of arm 0) are censored at day
# `control_censor_day` (181 = grace + 1); NULL follows them to the end of follow-up.

#' Surgery times of the hypothetical trial's surgery arm
#' @param s Planned surgery times of the rows (Inf if none); @param grace Grace period.
#' @param pool Observed surgery times within the grace period (donor dates).
#' @param seed Seed of the draw.
#' @return Numeric vector: own time if within grace, otherwise a donor date.
rct_surgery_times <- function(s, grace, pool, seed) {
  out <- s
  miss <- !(s <= grace)
  out[miss] <- withr::with_seed(seed, sample(pool, sum(miss), replace = TRUE),
                                .rng_kind = "L'Ecuyer-CMRG")
  out
}

#' Unweighted log HR of one hypothetical trial
#' @param dgm Fitted DGM; @param big Stacked base data; @param theta Treatment log-HR.
#' @param pool,grace,censor_day As above; @param seed Replicate seed.
#' @return Numeric log HR (arm 1 vs arm 0).
rct_estimate <- function(dgm, big, theta, pool, grace, censor_day, seed) {
  n <- nrow(big)
  rng <- draw_replicate_rng(seed, 2L * n)
  i1 <- seq_len(n)
  i0 <- n + i1
  d1 <- data.table::copy(big)
  d1[, s := rct_surgery_times(big$s, grace, pool, rng$seed_aux)]
  s1 <- simulate_outcomes(d1, dgm, theta, rng$u_event[i1], rng$u_cens[i1])
  d0 <- data.table::copy(big)[, s := Inf]
  s0 <- simulate_outcomes(d0, dgm, theta, rng$u_event[i0], rng$u_cens[i0])
  if (!is.null(censor_day)) {
    s0[, `:=`(status = as.integer(status == 1L & time <= censor_day), time = pmin(time, censor_day))]
  }
  dat <- data.table::rbindlist(list(s1[, .(time, status, arm = 1L)], s0[, .(time, status, arm = 0L)]))
  fit <- survival::coxph(survival::Surv(time, status) ~ arm, data = dat, ties = "efron")
  unname(stats::coef(fit)["arm"])
}

#' Large-sample value of the hypothetical-trial estimand at one theta
#' @param dgm Fitted DGM; @param data Base data; @param theta Treatment log-HR in the DGM.
#' @param pipeline Pipeline spec (grace, endpoint "hr"); @param cfg Full config (`rct`).
#' @param n_rep,n_target Large datasets to average and patients per dataset; @param seed Seed.
#' @return One-row data.table as `compute_truth()`: theta, truth, truth_mcse, n_truth, n_rep.
compute_truth_rct <- function(dgm, data, theta, pipeline, cfg, n_rep = cfg$truth$n_rep,
                              n_target = cfg$truth$n_target, seed = cfg$seed_truth) {
  checkmate::assert_number(theta)
  checkmate::assert_choice(endpoint_info(pipeline)$endpoint, "hr")
  big <- stack_copies(data, ceiling(n_target / nrow(data)))
  pool <- data$s[is.finite(data$s) & data$s <= pipeline$grace]
  cd <- cfg$rct$control_censor_day
  seeds <- make_crn_seeds(seed, n_rep)
  est <- vapply(seeds, function(sd) {
    rct_estimate(dgm, big, theta, pool, pipeline$grace, cd, sd)
  }, numeric(1))
  data.table::data.table(
    theta = theta, truth = mean(est),
    truth_mcse = if (n_rep > 1) stats::sd(est) / sqrt(n_rep) else NA_real_,
    n_truth = nrow(big), n_rep = n_rep
  )
}

#' Replace the CCW truth in nested outer replicates by the hypothetical-trial truth
#' Rebuilds each outer dataset from its seed, refits the prediction model (both are deterministic)
#' and recomputes the truth at the outer replicate's calibrated theta; the inner runs are reused.
#' @param nested_obs Outer replicates at theta* (a batch of `run_nested_batch()` output).
#' @param data Base data; @param dgm Fitted DGM; @param calibration From `calibrate_observed()`.
#' @param seeds Outer seeds; @param spec Dataset entry; @param cfg Full config.
#' @param level Nominal level.
#' @return `nested_obs` with truth_r, truth_r_mcse, bias_r, est_adj and the CI limits updated.
rct_nested_update <- function(nested_obs, data, dgm, calibration, seeds, spec, cfg, level = 0.95) {
  nc <- cfg$nested
  rows <- lapply(seq_len(nrow(nested_obs)), function(i) {
    r <- data.table::copy(nested_obs[i])
    if (is.na(r$est_adj)) return(r)
    rng <- draw_replicate_rng(seeds[r$rep], nrow(data))
    sim <- simulate_outcomes(data, dgm, calibration$theta_cal, rng$u_event, rng$u_cens)
    dgm_r <- fit_dgm(sim, spec$dgm, spec$censoring)
    sd <- nested_seeds(rng$seed_aux, nc$inner)
    tr <- tryCatch(
      compute_truth_rct(dgm_r, data, r$theta_r, spec$pipeline, cfg, n_rep = nc$truth_n_rep,
                        n_target = nc$truth_n_target, seed = sd$truth),
      error = function(e) NULL
    )
    if (is.null(tr)) {
      r[, c("truth_r", "truth_r_mcse", "bias_r", "est_adj") := NA_real_]
      return(r)
    }
    r[, `:=`(truth_r = tr$truth, truth_r_mcse = tr$truth_mcse)]
    r[, bias_r := mean_inner - truth_r]
    r[, est_adj := estimate - bias_r]
    r[, `:=`(ci_lo_np = est_adj + q_lo, ci_hi_np = est_adj + q_hi,
             ci_lo_np_refl = est_adj - q_hi, ci_hi_np_refl = est_adj - q_lo)]
    r[]
  })
  data.table::rbindlist(rows, fill = TRUE)
}

#' Hypothetical-trial truth at the scenario values (theta* and 0)
#' @param dgm Fitted DGM; @param data Base data; @param theta Scenario values.
#' @param pipeline Pipeline spec; @param cfg Full config.
#' @return data.table as `compute_truth_rct()`, one row per theta.
truth_rct_grid <- function(dgm, data, theta, pipeline, cfg) {
  data.table::rbindlist(lapply(theta, function(th) compute_truth_rct(dgm, data, th, pipeline, cfg)))
}
