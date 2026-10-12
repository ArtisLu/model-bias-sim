# Analysis pipelines under test: clone-censor-weight (CCW), one design per dataset.
#   grace_vs_none:    Maringe et al. (2020), surgery within grace vs no surgery within grace.
#   window_vs_window: Webster-Clark et al. (2025) tutorial, start within 90 vs within 30 days.

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
  long <- switch(spec$design,
    grace_vs_none = ccw_grace_vs_none(data, spec),
    window_vs_window = ccw_window_vs_window(data, spec)
  )
  fit_outcome(long, spec, robust)
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
