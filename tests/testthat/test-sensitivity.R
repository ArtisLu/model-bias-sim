# Small three-model chain (Weibull, exponential, generalized gamma), built once for these tests
sens_chain <- local({
  cfg_s <- config::get(file = test_path("..", "..", "config", "config.yml"), config = "sens_smoke")
  keys <- unlist(cfg_s$active_datasets)
  seeds <- make_crn_seeds(cfg_s$seed_master, 30)
  vseeds <- make_crn_seeds(cfg_s$seed_validation, 30)
  one <- function(ds) {
    sp <- cfg_s$datasets[[ds]]
    dgm <- fit_dgm(base_test, sp$dgm, sp$censoring)
    obs <- run_pipeline(base_test, sp$pipeline)
    th <- c(-0.7, 0)
    tt <- make_truth_table(rbindlist(lapply(th, function(t) {
      compute_truth(dgm, base_test, t, sp$pipeline, cfg_s, n_rep = 2, n_target = 1000)
    })), ds)
    run <- function(sd) {
      rbindlist(lapply(th, function(t) run_batch(base_test, dgm, t, sd, 1, sp$pipeline, 30)))[
        , dataset := ds]
    }
    sims <- run(seeds)
    adj <- adjust_simple(obs, sims, tt, ds)
    list(dgm = dgm, sims = sims, perf = summarise_performance(sims, tt), adj = adj,
         long = compare_performance(sims, tt, validate_adjustment(adj, run(vseeds), tt), NULL))
  }
  res <- stats::setNames(lapply(keys, one), keys)
  list(
    labels = vapply(cfg_s$datasets[keys], `[[`, character(1), "label"),
    dgms = stats::setNames(lapply(res, `[[`, "dgm"), paste0("dgm_", keys)),
    sims = rbindlist(lapply(res, `[[`, "sims")),
    perf = rbindlist(lapply(res, `[[`, "perf")),
    adj = rbindlist(lapply(res, `[[`, "adj")),
    long = rbindlist(stats::setNames(lapply(res, `[[`, "long"), paste0("perf_compare_", keys)),
                     idcol = "target")
  )
})

test_that("sensitivity profile: three models on the same data, nesting switched off", {
  cfg_s <- config::get(file = test_path("..", "..", "config", "config.yml"), config = "sens")
  expect_false(cfg_s$nested$enabled)
  expect_equal(unlist(cfg_s$active_datasets), c("maringe", "maringe_exp", "maringe_gg"))
  fam <- vapply(cfg_s$datasets[unlist(cfg_s$active_datasets)], function(x) x$dgm$family, "")
  expect_equal(unname(fam), c("weibull", "exponential", "gengamma"))
  # Same pipeline and covariates under every model
  expect_identical(cfg_s$datasets$maringe_exp$pipeline, cfg_s$datasets$maringe$pipeline)
  expect_identical(cfg_s$datasets$maringe_gg$pipeline, cfg_s$datasets$maringe$pipeline)
  # The default profile carries neither switch, so its stores stay valid
  expect_null(cfg_test$sensitivity)
  expect_null(cfg_test$nested$enabled)
})

test_that("model-fit table: nested models are ordered by likelihood and matched by dataset key", {
  ft <- sens_model_fit(sens_chain$dgms, sens_chain$labels)
  expect_equal(ft$model, c("Weibull", "Exponential", "Generalized gamma"))
  expect_equal(ft$n_par, ft$n_par[1] + c(0, -1, 1))
  expect_true(ft$loglik[2] <= ft$loglik[1] && ft$loglik[1] <= ft$loglik[3] + 1e-6)
  expect_true(is.na(ft$lr_p[1]) && all(ft$lr_p[2:3] >= 0 & ft$lr_p[2:3] <= 1))
  expect_equal(ft$aic, 2 * ft$n_par - 2 * ft$loglik)
})

test_that("comparison table reproduces the per-model adjustment and differences from the primary", {
  st <- sens_compare(sens_chain$adj, sens_chain$perf, sens_chain$long, sens_chain$sims, "maringe")
  expect_equal(st$dataset, c("maringe", "maringe_exp", "maringe_gg"))
  expect_equal(st$bias, sens_chain$adj$bias)
  expect_equal(st$d_bias, st$bias - st$bias[1])
  expect_equal(st$d_bias_mcse[1], 0)
  expect_true(all(st$d_bias_mcse[-1] > 0))
  expect_true(all(st$q_lo < 0 & st$q_hi > 0 & st$q_lo_mcse > 0))
  expect_true(all(is.finite(st$coverage_adj) & is.finite(st$type1_adj)))
  expect_error(sens_compare(sens_chain$adj, sens_chain$perf, sens_chain$long, sens_chain$sims, "x"))
})

test_that("quantile MCSE matches the spread of the quantile over repeated samples", {
  q <- withr::with_seed(1, replicate(400, stats::quantile(rnorm(2000), 0.025, type = 8)))
  m <- withr::with_seed(2, mean(replicate(50, quantile_mcse(rnorm(2000), 0.025))))
  expect_equal(m, stats::sd(q), tolerance = 0.15)
})

test_that("conservative combination: widest limits, estimate nearest to no effect", {
  ad <- data.table(
    dataset = c("a", "b"), est_adj = c(-0.5, -0.4), ci_lo_np = c(-1.0, -0.9),
    ci_hi_np = c(0.0, 0.1), null_q_lo = c(-0.52, -0.50), null_q_hi = c(0.55, 0.57),
    z_crit_par = c(2.01, 2.04), p_np = c(0.04, 0.06), reject_np = c(TRUE, FALSE)
  )
  cs <- sens_conservative(ad)
  expect_equal(cs$est_adj, -0.4)
  expect_equal(c(cs$ci_lo_np, cs$ci_hi_np), c(-1.0, 0.1))
  expect_equal(c(cs$null_q_lo, cs$null_q_hi), c(-0.52, 0.57))
  expect_equal(cs$z_crit_par, 2.04)
  expect_equal(cs$p_np, 0.06)
  expect_false(cs$reject_np)
  expect_equal(sens_conservative(ad, keep = "a")$est_adj, -0.5)
})

test_that("sensitivity figures build", {
  expect_s3_class(plot_sens_fit(base_test, sens_chain$dgms, sens_chain$labels), "ggplot")
  expect_s3_class(plot_sens_distribution(sens_chain$sims, sens_chain$perf, sens_chain$labels),
                  "ggplot")
})
