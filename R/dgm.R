# Reference data-generating model: PH model for the outcome with treatment as a time-varying
# exposure and a Weibull (default) or 3-parameter generalized gamma baseline hazard (sensitivity
# analysis, paper Section 4.6), plus a Weibull PH model for loss to follow-up fitted with the
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

#' Baseline of the generalized gamma (flexsurv parametrisation, mu = 0): cumulative hazard
#' H0(t) = -log S0(t), log hazard, and the inverse of H0. With Q = 1 it is the Weibull baseline
#' H0(t) = t^(1 / sigma); Q -> 0 is the log-normal. The scale is carried by the PH intercept.
#' @param t Times (may be Inf); @param sigma Scale > 0; @param Q Shape (any sign).
#' @return Numeric vector.
gg_cumhaz0 <- function(t, sigma, Q) {
  w <- log(t) / sigma
  if (abs(Q) < 1e-6) return(-stats::pnorm(w, lower.tail = FALSE, log.p = TRUE))
  k <- Q^-2
  x <- k * exp(Q * w)
  -stats::pgamma(x, k, lower.tail = Q < 0, log.p = TRUE)
}

#' @rdname gg_cumhaz0
gg_loghaz0 <- function(t, sigma, Q) {
  w <- log(t) / sigma
  if (abs(Q) < 1e-6) {
    return(stats::dnorm(w, log = TRUE) - log(sigma * t) -
             stats::pnorm(w, lower.tail = FALSE, log.p = TRUE))
  }
  k <- Q^-2
  x <- k * exp(Q * w)
  logf <- log(abs(Q)) - log(sigma * t) + k * log(k) - lgamma(k) + k * (Q * w - exp(Q * w))
  logf - stats::pgamma(x, k, lower.tail = Q < 0, log.p = TRUE)
}

#' @rdname gg_cumhaz0
#' @param H Baseline cumulative hazard values.
gg_time0 <- function(H, sigma, Q) {
  if (abs(Q) < 1e-6) return(exp(sigma * stats::qnorm(-H, lower.tail = FALSE, log.p = TRUE)))
  k <- Q^-2
  x <- stats::qgamma(-H, k, lower.tail = Q < 0, log.p = TRUE)
  exp(sigma * log(Q^2 * x) / Q)
}

#' Negative log-likelihood of the PH generalized gamma model with time-varying treatment
#' Time is divided by `tau` (a fixed time unit, here the length of follow-up) so that the PH
#' intercept stays moderate; the baseline is H0(t / tau).
#' @param par c(alpha, log_sigma, Q, gamma, theta); @param X,time,status,s As `dgm_negloglik()`.
#' @param tau Time unit.
#' @return Scalar (Inf if the likelihood is not finite).
gg_negloglik <- function(par, X, time, status, s, tau) {
  p <- ncol(X)
  sigma <- exp(par[2])
  Q <- par[3]
  eta <- par[1] + drop(X %*% par[3 + seq_len(p)])
  theta <- par[p + 4]
  on_trt <- s < time
  H0t <- gg_cumhaz0(time / tau, sigma, Q)
  H0s <- gg_cumhaz0(pmin(s, time) / tau, sigma, Q)
  cumhaz <- exp(eta) * (H0s + on_trt * exp(theta) * (H0t - H0s))
  loghaz <- gg_loghaz0(time / tau, sigma, Q) - log(tau) + eta + theta * on_trt
  nll <- -sum(status * loghaz - cumhaz)
  if (is.finite(nll)) nll else Inf
}

#' Fit the PH generalized gamma model with time-varying treatment by maximum likelihood
#' Started at the Weibull fit (the Weibull is the Q = 1 member), so the likelihood cannot be lower.
#' @param data Base data; @param covariates Covariate names; @param status Event indicator.
#' @return List like `fit_weibull_tv()` with family "gengamma", sigma, Q (rho is NA) and the
#'   Weibull log-likelihood `loglik_weibull` for a likelihood-ratio check.
fit_gengamma_tv <- function(data, covariates, status) {
  wb <- fit_weibull_tv(data, covariates, status)
  X <- dgm_design(data, covariates)
  p <- ncol(X)
  tau <- max(data$admin_time)
  # Weibull t^rho = tau^rho (t / tau)^rho: its intercept moves by rho * log(tau)
  start <- c(wb$alpha + wb$rho * log(tau), -log(wb$rho), 1, wb$gamma, wb$theta_hat)
  fit <- stats::optim(
    start, gg_negloglik, X = X, time = data$time, status = status, s = data$s, tau = tau,
    method = "BFGS", hessian = TRUE, control = list(maxit = 2000, reltol = 1e-12)
  )
  if (fit$convergence != 0) cli::cli_abort("Generalized gamma fit did not converge.")
  vc <- tryCatch(solve(fit$hessian), error = function(e) matrix(NA_real_, p + 4, p + 4))
  nm <- c("alpha", "log_sigma", "Q", colnames(X), "theta")
  dimnames(vc) <- list(nm, nm)
  list(
    family = "gengamma",
    alpha = unname(fit$par[1]), sigma = unname(exp(fit$par[2])), Q = unname(fit$par[3]), tau = tau,
    rho = NA_real_,
    gamma = stats::setNames(unname(fit$par[3 + seq_len(p)]), colnames(X)),
    theta_hat = unname(fit$par[p + 4]),
    se = stats::setNames(sqrt(pmax(diag(vc), 0)), nm),
    covariates = covariates,
    loglik = -fit$value, loglik_weibull = wb$loglik,
    n_events = sum(status)
  )
}

#' Fit the reference DGM (outcome and censoring models) to the base data
#' @param data Base data from `prepare_base_data()`.
#' @param spec `cfg$dgm` (`family`: weibull, default, or gengamma);
#' @param cens_spec `cfg$censoring` (model: none, fixed or weibull).
#' @return Outcome-model list (see `fit_weibull_tv()`, `fit_gengamma_tv()`) with element `cens`.
fit_dgm <- function(data, spec, cens_spec = list(model = "fixed")) {
  checkmate::assert_data_table(data)
  checkmate::assert_choice(cens_spec$model, c("none", "fixed", "weibull"))
  family <- if (is.null(spec$family)) "weibull" else spec$family
  checkmate::assert_choice(family, c("weibull", "gengamma"))
  dgm <- if (family == "gengamma") {
    fit_gengamma_tv(data, spec$covariates, data$status)
  } else {
    c(list(family = "weibull"), fit_weibull_tv(data, spec$covariates, data$status))
  }
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

#' Invert the outcome model's cumulative hazard (any family) with treatment switching at time s
#' @param dgm Fitted DGM; @param eta,theta,s,u As `invert_weibull_tv()`.
#' @return Simulated times.
invert_dgm_tv <- function(dgm, eta, theta, s, u) {
  if (is.null(dgm$family) || dgm$family == "weibull") {
    return(invert_weibull_tv(eta, dgm$rho, theta, s, u))
  }
  E <- -log(u)
  H0s <- gg_cumhaz0(s / dgm$tau, dgm$sigma, dgm$Q)
  h_pre <- exp(eta) * H0s
  H0t <- data.table::fifelse(E <= h_pre, E / exp(eta), H0s + (E - h_pre) / exp(eta + theta))
  dgm$tau * gg_time0(H0t, dgm$sigma, dgm$Q)
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
  t_ev <- invert_dgm_tv(dgm, dgm$alpha + drop(X %*% dgm$gamma), theta, data$s, u)
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
