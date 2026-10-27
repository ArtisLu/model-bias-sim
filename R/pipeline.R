# Analysis pipelines under test: clone-censor-weight (CCW), one design per dataset.
#   grace_vs_none:    Maringe et al. (2020), surgery within grace vs no surgery within grace.
#   window_vs_window: Webster-Clark et al. (2025) tutorial, start within 90 vs within 30 days.
# Endpoint (spec$endpoint): "hr" = log hazard ratio of the weighted Cox model (default);
# "surv12" = difference in survival at spec$surv_time, from IPC-weighted proportions (Maringe only).

#' Effect scale of an endpoint: labels and the direction that favours arm 1
#' @param spec Pipeline spec.
#' @return list(endpoint, scale, short, favours_low, null_value). `favours_low` is TRUE when a
#'   smaller estimate favours arm 1 (log HR); FALSE for a survival difference.
endpoint_info <- function(spec) {
  ep <- if (is.null(spec$endpoint)) "hr" else spec$endpoint
  checkmate::assert_choice(ep, c("hr", "surv12"))
  switch(ep,
    hr = list(endpoint = ep, scale = "log hazard ratio", short = "log HR", favours_low = TRUE),
    surv12 = list(endpoint = ep, scale = "difference in 12-month survival",
                  short = "survival difference", favours_low = FALSE)
  )
}

#' Split clones into intervals at given breaks
#' @param clones Clone rows of one arm.
#' @param breaks Interval start times (first = 0); last interval is open.
#' @return Long data.table with k, start, stop, ev (event in interval).
split_clones <- function(clones, breaks) {
  stops <- c(breaks[-1], Inf)
  n_k <- findInterval(clones$end, breaks, left.open = TRUE)
  long <- clones[rep(seq_len(.N), n_k)]
  long[, k := sequence(n_k)]
  long[, `:=`(start = breaks[k], stop = pmin(stops[k], end))]
  long[, last := k == rep(n_k, n_k)]
  long[, ev := as.integer(last & event == 1L)]
  long[]
}

#' Weighted Cox model for the strategy contrast
#' @param long Stacked long data (id, arm, start, stop, ev, w); @param spec Pipeline spec.
#' @param robust Robust (cluster) SE; FALSE only for large-sample truth runs.
#' @return One-row data.table in the pipeline contract format.
fit_outcome <- function(long, spec, robust) {
  fit <- if (robust) {
    survival::coxph(
      survival::Surv(start, stop, ev) ~ arm,
      data = long, weights = w, cluster = id, ties = "efron"
    )
  } else {
    survival::coxph(survival::Surv(start, stop, ev) ~ arm, data = long, weights = w, ties = "efron")
  }
  est <- unname(stats::coef(fit)["arm"])
  se <- unname(sqrt(stats::vcov(fit)["arm", "arm"]))
  z <- stats::qnorm(1 - (1 - spec$conf_level) / 2)
  data.table::data.table(
    estimate = est, se = se, ci_lo = est - z * se, ci_hi = est + z * se,
    p_value = 2 * stats::pnorm(-abs(est / se)), n_events = sum(long$ev),
    status = "ok", message = NA_character_
  )
}

#' Run the CCW pipeline (same function for real and simulated data)
#' @param data Patient-level data in the common format; @param spec Pipeline spec.
#' @param robust Robust (cluster) SE; FALSE only for large-sample truth runs.
#' @return One-row data.table: estimate (log-HR, arm 1 vs arm 0), se, ci_lo, ci_hi,
#'   p_value, n_events, status, message.
run_pipeline <- function(data, spec, robust = TRUE) {
  checkmate::assert_data_table(data)
  checkmate::assert_choice(spec$design, c("grace_vs_none", "window_vs_window"))
  ep <- endpoint_info(spec)$endpoint
  long <- switch(spec$design,
    grace_vs_none = ccw_grace_vs_none(data, spec),
    window_vs_window = ccw_window_vs_window(data, spec)
  )
  switch(ep,
    hr = fit_outcome(long, spec, robust),
    surv12 = fit_survival_difference(long, spec)
  )
}

#' Difference in survival at a fixed time from the CCW clones (binary endpoint)
#' Hajek estimator per arm: sum(w * alive) / sum(w) over clones still uncensored at the time, with
#' the weight each clone has at the end of its follow-up (the same IPC weights as the Cox model).
#' The SE sums the arm-wise influence functions within patient (the two clones of a patient are
#' dependent), treating the weights as known, like the cluster-robust SE of the Cox model.
#' @param long Stacked long data (id, arm, start, stop, ev, w); @param spec Pipeline spec.
#' @return One-row data.table in the pipeline contract format; estimate = S(arm 1) - S(arm 0).
fit_survival_difference <- function(long, spec) {
  checkmate::assert_number(spec$surv_time, lower = 0)
  checkmate::assert_choice(spec$design, "grace_vs_none")
  t_s <- spec$surv_time
  # One row per clone: its last interval carries the end of follow-up and the final weight
  last <- long[long[, .I[.N], by = .(id, arm)]$V1]
  last[, `:=`(
    observed = as.integer(ev == 1L | stop >= t_s),
    alive = as.integer(stop >= t_s)
  )]
  last <- last[observed == 1L]
  arm_fit <- function(a) {
    x <- last[arm == a]
    s <- sum(x$w * x$alive) / sum(x$w)
    x[, infl := w * (alive - s) / sum(w)]
    list(s = s, infl = x[, .(id, infl)])
  }
  f1 <- arm_fit(1L)
  f0 <- arm_fit(0L)
  infl <- merge(f1$infl, f0$infl, by = "id", all = TRUE, suffixes = c("1", "0"))
  infl[is.na(infl1), infl1 := 0][is.na(infl0), infl0 := 0]
  est <- f1$s - f0$s
  se <- sqrt(sum((infl$infl1 - infl$infl0)^2))
  z <- stats::qnorm(1 - (1 - spec$conf_level) / 2)
  data.table::data.table(
    estimate = est, se = se, ci_lo = est - z * se, ci_hi = est + z * se,
    p_value = 2 * stats::pnorm(-abs(est / se)), n_events = sum(long$ev),
    status = "ok", message = NA_character_
  )
}

# ---- Maringe design ------------------------------------------------------------

#' Clone into "surgery within grace" (arm 1) and "no surgery within grace" (arm 0)
#' @param data Patient-level data; @param grace Grace period (days).
#' @return Clone rows: arm, end, event, art_cens.
build_clones <- function(data, grace) {
  ctrl <- data.table::copy(data)[, arm := 0L]
  ctrl[, art_cens := as.integer(trt_time <= grace)]
  ctrl[, `:=`(
    end = data.table::fifelse(art_cens == 1L, trt_time, time),
    event = data.table::fifelse(art_cens == 1L, 0L, status)
  )]
  surg <- data.table::copy(data)[, arm := 1L]
  # Deaths before grace without surgery stay in the surgery arm (Maringe et al. 2020)
  surg[, art_cens := as.integer(trt_time > grace & time > grace)]
  surg[, `:=`(
    end = data.table::fifelse(art_cens == 1L, grace, time),
    event = data.table::fifelse(art_cens == 1L, 0L, status)
  )]
  data.table::rbindlist(list(ctrl, surg))
}

#' Inverse probability of censoring weights for the control arm (Cox censoring model)
#' @param clones Control-arm clones; @param spec Pipeline spec.
#' @return Long data split at censoring times, with weight `w` = 1 / P(uncensored at start).
weight_control <- function(clones, spec) {
  f <- stats::reformulate(spec$cens_covariates, response = "survival::Surv(end, art_cens)")
  fit <- survival::coxph(f, data = clones, ties = "breslow")
  bh <- survival::basehaz(fit, centered = FALSE)
  clones[, lp := stats::predict(fit, type = "lp", reference = "zero")]
  # Split at the observed censoring (surgery) times: weights change only there
  ctimes <- sort(unique(clones[art_cens == 1L, end]))
  brk <- unique(c(0, ctimes[ctimes < spec$grace], spec$grace))
  long <- split_clones(clones, brk)
  H0 <- stats::stepfun(bh$time, c(0, bh$hazard))
  # Right-continuous: those censored at a break have left the risk set
  long[, w := exp(H0(pmin(start, spec$grace)) * exp(lp))]
  long[]
}

#' Inverse probability of censoring weights for the surgery arm
#' @param clones Surgery-arm clones; @param data Patient data; @param spec Pipeline spec.
#' @return Long data split at grace, with weight `w`.
weight_surgery <- function(clones, data, spec) {
  long <- split_clones(clones, c(0, spec$grace))
  # Censoring happens only at grace end, for patients alive and not yet operated
  alive <- data[time > spec$grace]
  alive[, op := as.integer(trt_time <= spec$grace)]
  f <- stats::reformulate(spec$cens_covariates, response = "op")
  fit <- stats::glm(f, family = stats::binomial(), data = alive)
  alive[, p_op := stats::predict(fit, type = "response")]
  long[, w := 1]
  long[alive, on = "id", p_op := i.p_op]
  long[start >= spec$grace, w := 1 / p_op]
  long[]
}

#' CCW long data for the Maringe design
#' @param data Patient-level data; @param spec Pipeline spec.
#' @return Stacked long data (id, arm, start, stop, ev, w).
ccw_grace_vs_none <- function(data, spec) {
  clones <- build_clones(data, spec$grace)
  long0 <- weight_control(clones[arm == 0L], spec)
  long1 <- weight_surgery(clones[arm == 1L], data, spec)
  cols <- c("id", "arm", "start", "stop", "ev", "w")
  data.table::rbindlist(list(long0[, ..cols], long1[, ..cols]))
}

# ---- Webster-Clark tutorial design ----------------------------------------------

#' One "start treatment within `win` days" clone set, censored and weighted as in the tutorial
#' @param data Patient-level data; @param win Window (days); @param spec Pipeline spec.
#' @return Long data (id, start, stop, ev, w).
ccw_window_arm <- function(data, win, spec) {
  cl <- data.table::copy(data)
  # Tutorial 01: censored at win if followed beyond win without having started
  cl[, nostart := as.integer(time > win & trt_time > win)]
  cl[, `:=`(
    end = data.table::fifelse(nostart == 1L, win, time),
    event = data.table::fifelse(nostart == 1L, 0L, status)
  )]
  long <- split_clones(cl, c(0, win))
  # Tutorial 02: weight 1 in [0, win]; after win, up-weight only recent starters
  after <- long[start >= win]
  model_rows <- rbind(
    cl[nostart == 1L],
    cl[nostart == 0L & time > win & trt_time >= win - spec$recent_days]
  )
  model_rows[, uncens := 1L - nostart]
  f <- stats::reformulate(spec$cens_covariates, response = "uncens")
  fit <- stats::glm(f, family = stats::binomial(), data = model_rows)
  after[, p_unc := stats::predict(fit, newdata = after, type = "response")]
  after[, w := data.table::fifelse(trt_time >= win - spec$recent_days, 1 / p_unc, 1)]
  long[, w := 1]
  long[after, on = c("id", "k"), w := i.w]
  long[]
}

#' CCW long data for the tutorial design (arm 1 = longer window vs arm 0 = shorter window)
#' @param data Patient-level data; @param spec Pipeline spec.
#' @return Stacked long data (id, arm, start, stop, ev, w).
ccw_window_vs_window <- function(data, spec) {
  win <- unlist(spec$windows)
  cols <- c("id", "start", "stop", "ev", "w")
  long0 <- ccw_window_arm(data, win[1], spec)[, ..cols][, arm := 0L]
  long1 <- ccw_window_arm(data, win[2], spec)[, ..cols][, arm := 1L]
  data.table::rbindlist(list(long0, long1), use.names = TRUE)
}
