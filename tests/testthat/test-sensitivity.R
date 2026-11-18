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
    # Same truth seeds under every model, as in the pipeline
    tseeds <- make_crn_seeds(cfg_s$seed_truth, 3)
    rows <- rbindlist(lapply(th, function(t) {
      rbindlist(lapply(1:3, function(j) {
        truth_point(dgm, base_test, t, j, tseeds, sp$pipeline, cfg_s, n_target = 1000)
      }))
    }))
    tt <- make_truth_table(summarise_truth(rows), ds)
    run <- function(sd) {
      rbindlist(lapply(th, function(t) run_batch(base_test, dgm, t, sd, 1, sp$pipeline, 30)))[
        , dataset := ds]
    }
    sims <- run(seeds)
    adj <- adjust_simple(obs, sims, tt, ds)
    list(dgm = dgm, sims = sims, perf = summarise_performance(sims, tt), adj = adj,
         rows = rows[, dataset := ds][],
         long = compare_performance(sims, tt, validate_adjustment(adj, run(vseeds), tt), NULL))
  }
  res <- stats::setNames(lapply(keys, one), keys)
  list(
    labels = vapply(cfg_s$datasets[keys], `[[`, character(1), "label"),
    dgms = stats::setNames(lapply(res, `[[`, "dgm"), paste0("dgm_", keys)),
    sims = rbindlist(lapply(res, `[[`, "sims")),
    perf = rbindlist(lapply(res, `[[`, "perf")),
    adj = rbindlist(lapply(res, `[[`, "adj")),
    rows = rbindlist(lapply(res, `[[`, "rows")),
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
  expect_false(any(st$d_paired))
})

test_that("pairing the true values by seed gives a valid and smaller MCSE of the bias difference", {
  ind <- sens_compare(sens_chain$adj, sens_chain$perf, sens_chain$long, sens_chain$sims, "maringe")
  prd <- sens_compare(sens_chain$adj, sens_chain$perf, sens_chain$long, sens_chain$sims, "maringe",
                      truth_rows = sens_chain$rows)
  expect_equal(prd$d_bias, ind$d_bias)
  expect_equal(prd$d_paired, c(FALSE, TRUE, TRUE))
  expect_true(all(prd$d_bias_mcse[-1] > 0 & prd$d_bias_mcse[-1] < ind$d_bias_mcse[-1]))
  # Direct computation for one model
  r <- sens_chain$rows
  th <- sens_chain$perf[scenario == "observed", setNames(theta, dataset)]
  a <- r[dataset == "maringe_gg" & abs(theta - th[["maringe_gg"]]) < 1e-10][order(seed_id)]
  b <- r[dataset == "maringe" & abs(theta - th[["maringe"]]) < 1e-10][order(seed_id)]
  expect_equal(paired_truth_mcse(r, "maringe_gg", th[["maringe_gg"]], "maringe", th[["maringe"]]),
               sd(a$estimate - b$estimate) / sqrt(3))
  expect_true(is.na(paired_truth_mcse(NULL, "maringe_gg", 0, "maringe", 0)))
})

test_that("limit differences: paired quantile differences with bootstrap MCSE", {
  ld <- sens_limit_differences(sens_chain$sims, sens_chain$perf, "maringe", n_boot = 50)
  expect_equal(nrow(ld), 8)
  expect_setequal(ld$limit, c("ci_lo", "ci_hi", "test_lo", "test_hi"))
  expect_true(all(ld$diff_mcse > 0))
  # Against the primary model itself the difference and its MCSE are zero
  self <- copy(sens_chain$sims[dataset == "maringe"])[, dataset := "copy"]
  pf <- rbind(sens_chain$perf, copy(sens_chain$perf[dataset == "maringe"])[, dataset := "copy"])
  z <- sens_limit_differences(rbind(sens_chain$sims, self), pf, "maringe", n_boot = 20)
  expect_equal(z[dataset == "copy", diff], rep(0, 4))
  expect_equal(z[dataset == "copy", diff_mcse], rep(0, 4))
})

test_that("verdict: negligible, unclear or relevant against the threshold, allowing for MCSE", {
  st <- data.table(dataset = c("p", "m"), d_bias = c(0, 0.03), d_bias_mcse = c(0, 0.005),
                   emp_se = c(0.3, 0.3))
  ld <- data.table(dataset = "m", limit = c("ci_lo", "ci_hi", "test_lo"),
                   diff = c(0.004, 0.016, -0.001), diff_mcse = c(0.001, 0.004, 0.008))
  v <- sens_verdict(st, ld)
  expect_equal(v$threshold[1], 0.015)
  expect_equal(v$detectable, c(TRUE, TRUE, TRUE, FALSE))
  # 0.03 - 0.01 > 0.015; 0.004 + 0.002 < 0.015; 0.016 within 0.015 +/- 0.008; 0.001 + 0.016 > 0.015
  expect_equal(v$status, c("relevant", "negligible", "unclear", "unclear"))
})
