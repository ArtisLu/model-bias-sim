# Reference data-generating model: Weibull PH for the outcome with treatment as a
# time-varying exposure, plus a Weibull PH model for loss to follow-up fitted with the
# event/censoring indicator switched (paper Section 2.3). Simulation by inversion (CRN).

#' Design matrix of DGM covariates (no intercept)
#' @param data Base data.
#' @param covariates Covariate names.
#' @return Numeric matrix.
dgm_design <- function(data, covariates) {
  f <- stats::reformulate(covariates)
  stats::model.matrix(f, data = data)[, -1, drop = FALSE]
}

#' Negative log-likelihood of the Weibull PH model with time-varying treatment
#' @param par c(alpha, log_rho, gamma, theta).
#' @param X Design matrix; @param time,status,s Follow-up, event, planned treatment time.
#' @return Scalar.
dgm_negloglik <- function(par, X, time, status, s) {
  p <- ncol(X)
  rho <- exp(par[2])
  eta <- par[1] + drop(X %*% par[2 + seq_len(p)])
  theta <- par[p + 3]
  on_trt <- s < time
  cumhaz <- exp(eta) *
    (pmin(s, time)^rho + on_trt * exp(theta) * (time^rho - pmin(s, time)^rho))
  loghaz <- log(rho) + (rho - 1) * log(time) + eta + theta * on_trt
  -sum(status * loghaz - cumhaz)
}

#' Fit a Weibull PH model with time-varying treatment by maximum likelihood
#' @param data Base data; @param covariates Covariate names.
#' @param status Event indicator (outcome, or loss to follow-up for the censoring model).
#' @return List: alpha, rho, gamma, theta_hat, se, covariates, loglik, n_events.
fit_weibull_tv <- function(data, covariates, status) {
  X <- dgm_design(data, covariates)
  start <- c(log(sum(status) / sum(data$time)), 0, rep(0, ncol(X)), 0)
  fit <- stats::optim(
    start, dgm_negloglik,
    X = X, time = data$time, status = status, s = data$s,
    method = "BFGS", hessian = TRUE, control = list(maxit = 1000, reltol = 1e-12)
  )
  if (fit$convergence != 0) cli::cli_abort("Weibull fit did not converge.")
  p <- ncol(X)
  vc <- solve(fit$hessian)
  nm <- c("alpha", "log_rho", colnames(X), "theta")
  dimnames(vc) <- list(nm, nm)
  list(
    alpha = fit$par[1],
    rho = exp(fit$par[2]),
    gamma = stats::setNames(fit$par[2 + seq_len(p)], colnames(X)),
    theta_hat = fit$par[p + 3],
    se = stats::setNames(sqrt(diag(vc)), nm),
    covariates = covariates,
    loglik = -fit$value,
    n_events = sum(status)
  )
}

#' Fit the reference DGM (outcome and censoring models) to the base data
#' @param data Base data from `prepare_base_data()`.
#' @param spec `cfg$dgm`; @param cens_spec `cfg$censoring` (model: none, fixed or weibull).
#' @return Outcome-model list (see `fit_weibull_tv()`) with element `cens`.
fit_dgm <- function(data, spec, cens_spec = list(model = "fixed")) {
  checkmate::assert_data_table(data)
  checkmate::assert_choice(cens_spec$model, c("none", "fixed", "weibull"))
  dgm <- fit_weibull_tv(data, spec$covariates, data$status)
  # Censoring model: loss to follow-up is the event; outcome events and the horizon censor it
  dgm$cens <- if (cens_spec$model == "weibull") {
    c(list(model = "weibull"), fit_weibull_tv(data, cens_spec$covariates, data$cens_event))
  } else {
    list(model = cens_spec$model)
  }
  dgm
}

#' Invert the piecewise Weibull cumulative hazard (switch to exp(theta) at time s)
#' @param eta Linear predictor incl. intercept; @param rho Shape; @param theta Log-HR on treatment.
#' @param s Planned treatment time; @param u Uniforms.
#' @return Simulated times; the same u gives times monotone in theta.
invert_weibull_tv <- function(eta, rho, theta, s, u) {
  E <- -log(u)
  h_pre <- exp(eta) * s^rho
  data.table::fifelse(
    E <= h_pre,
    (E / exp(eta))^(1 / rho),
    (s^rho + (E - h_pre) / exp(eta + theta))^(1 / rho)
  )
}

#' Simulate event (and censoring) times for fixed baseline data and treatment plans
#' @param data Base data (covariates and planned treatment time `s` are kept).
#' @param dgm Fitted DGM from `fit_dgm()`.
#' @param theta Treatment log-HR used for simulation.
#' @param u,u_cens Uniforms for event and censoring times, one per patient (CRN).
#' @return Copy of `data` with simulated `time`, `status`, `cens_time`, `trt_time`, `cens_event`.
simulate_outcomes <- function(data, dgm, theta, u, u_cens = NULL) {
  checkmate::assert_numeric(u, len = nrow(data), lower = 0, upper = 1)
  X <- dgm_design(data, dgm$covariates)
  t_ev <- invert_weibull_tv(dgm$alpha + drop(X %*% dgm$gamma), dgm$rho, theta, data$s, u)
  out <- data.table::copy(data)
  cm <- if (is.null(dgm$cens)) "fixed" else dgm$cens$model
  if (cm == "weibull") {
    checkmate::assert_numeric(u_cens, len = nrow(data), lower = 0, upper = 1)
    cs <- dgm$cens
    Xc <- dgm_design(data, cs$covariates)
    # Censoring hazard also changes at treatment start; theta_c is kept as fitted
    eta_c <- cs$alpha + drop(Xc %*% cs$gamma)
    t_c <- invert_weibull_tv(eta_c, cs$rho, cs$theta_hat, data$s, u_cens)
    out[, cens_time := pmin(t_c, admin_time)]
  }
  out[, `:=`(time = pmin(t_ev, cens_time), status = as.integer(t_ev <= cens_time))]
  out[, trt_time := data.table::fifelse(s <= time, s, Inf)]
  # Simulated loss to follow-up, so a refit of the censoring model (nested simulation) sees it
  out[, cens_event := as.integer(status == 0L & time < admin_time)]
  out[]
}
