# Reference data-generating model: Weibull PH for the outcome with treatment as a
# time-varying exposure; outcome simulation by inversion (CRN-compatible).

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

#' Fit the reference DGM to the base data
#' @param data Base data from `prepare_base_data()`.
#' @param spec `cfg$dgm`.
#' @return List: alpha, rho, gamma, theta_hat, vcov, covariates, convergence.
fit_dgm <- function(data, spec) {
  checkmate::assert_data_table(data)
  X <- dgm_design(data, spec$covariates)
  start <- c(log(sum(data$status) / sum(data$time)), 0, rep(0, ncol(X)), 0)
  fit <- stats::optim(
    start, dgm_negloglik,
    X = X, time = data$time, status = data$status, s = data$s,
    method = "BFGS", hessian = TRUE, control = list(maxit = 1000, reltol = 1e-12)
  )
  if (fit$convergence != 0) cli::cli_abort("DGM fit did not converge.")
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
    covariates = spec$covariates,
    loglik = -fit$value
  )
}

#' Simulate event times for fixed baseline data and treatment plans
#' @param data Base data (covariates, planned treatment time `s` and `cens_time` are kept).
#' @param dgm Fitted DGM from `fit_dgm()`.
#' @param theta Treatment log-HR used for simulation.
#' @param u Uniforms, one per patient (common random numbers).
#' @return Copy of `data` with simulated `time`, `status`, `trt_time`.
simulate_outcomes <- function(data, dgm, theta, u) {
  checkmate::assert_numeric(u, len = nrow(data), lower = 0, upper = 1)
  X <- dgm_design(data, dgm$covariates)
  eta <- dgm$alpha + drop(X %*% dgm$gamma)
  rho <- dgm$rho
  E <- -log(u)
  # Inversion of the piecewise cumulative hazard: same u gives monotone T in theta
  h_pre <- exp(eta) * data$s^rho
  t_ev <- data.table::fifelse(
    E <= h_pre,
    (E / exp(eta))^(1 / rho),
    (data$s^rho + (E - h_pre) / exp(eta + theta))^(1 / rho)
  )
  out <- data.table::copy(data)
  out[, `:=`(time = pmin(t_ev, cens_time), status = as.integer(t_ev <= cens_time))]
  out[, trt_time := data.table::fifelse(s <= time, s, Inf)]
  out[]
}
