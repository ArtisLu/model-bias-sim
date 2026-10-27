# _targets.R - modelling-bias QBA simulation (CCW), paper v0.65: inverse calibration, simple and
# nested non-parametric adjustment, for the analyses in cfg$active_datasets (default: Maringe, log
# HR and 12-month survival difference), pilot study for the outer / inner split (Appendix) and the
# hypothetical-trial target for a non-consistent estimator (Section 4.6, HR analysis only).
# Run with targets::tar_make(); smoke run: Sys.setenv(R_CONFIG_ACTIVE = "smoke");
# two-dataset comparison: "two_datasets"; generalized gamma sensitivity analysis: "gengamma".
library(targets)
library(tarchetypes)

tar_source("R")
cfg <- config::get(file = "config/config.yml")
# config attaches the absolute path of the file; dropping it keeps the store valid when the
# project is moved or cloned elsewhere (every target depends on cfg)
attributes(cfg) <- list(names = names(cfg))
# Optional CI comparison (bootstrap-t, confidence belt) in its own file: toggling it keeps cfg
cfg_belt <- config::get(file = "config/belt.yml")
attributes(cfg_belt) <- list(names = names(cfg_belt))
extras <- isTRUE(cfg_belt$enabled)
active <- unlist(cfg$active_datasets)
# Hazard-ratio analyses only: the comparison figures use the HR scale
is_hr <- vapply(active, function(d) {
  endpoint_info(cfg$datasets[[d]]$pipeline)$endpoint == "hr"
}, logical(1))
ds_labels <- vapply(cfg$datasets[active[is_hr]], `[[`, character(1), "label")
# Outputs of the default (paper) profile go to reports/; other profiles get their own folder so a
# precision run cannot overwrite the paper figures
profile <- Sys.getenv("R_CONFIG_ACTIVE", "default")
out_dir <- if (profile == "default") "reports" else file.path("reports", profile)
# Report on the Maringe case study (into out_dir): not in smoke runs, only if Quarto is installed.
# Kept out of cfg on purpose: every target depends on cfg, so a new key would rerun everything.
do_report <- "maringe" %in% active && profile != "smoke" &&
  requireNamespace("quarto", quietly = TRUE) &&
  isTRUE(tryCatch(nzchar(quarto::quarto_path()), error = function(e) FALSE))
# Parallel workers need crew; without it the pipeline runs sequentially. QBA_WORKERS overrides
# cfg$workers without editing the config (every target depends on cfg, so editing it reruns all).
n_workers <- as.integer(Sys.getenv("QBA_WORKERS", cfg$workers))
use_crew <- n_workers > 1 && requireNamespace("crew", quietly = TRUE)

tar_option_set(
  packages = c("data.table", "survival", "checkmate"),
  seed = cfg$seed_master,
  controller = if (use_crew) crew::crew_controller_local(workers = n_workers) else NULL
)

# Optional CI comparison per dataset (config/belt.yml, enabled: true)
extra_targets <- if (extras) {
  list(
    tar_target(belt_theta, belt_theta_grid(cfg_belt$datasets[[ds]]$theta_extra, theta_sim)),
    tar_target(
      belt_truth,
      compute_truth(dgm, base_data, belt_theta, cfg$datasets[[ds]]$pipeline, cfg,
                    n_rep = cfg_belt$truth_n_rep),
      pattern = map(belt_theta)
    ),
    tar_target(
      belt_sims,
      run_batch(base_data, dgm, belt_theta, crn_seeds[seq_len(min(cfg_belt$B, cfg$B))],
                belt_batch_id, cfg$datasets[[ds]]$pipeline, cfg$batch_size)[, dataset := ds],
      pattern = cross(belt_theta, belt_batch_id)
    ),
    # Column names as strings: tar_map() renames every symbol that matches a target (truth)
    tar_target(belt, fit_belt(
      rbind(sim_results, belt_sims, fill = TRUE),
      rbind(truth_table[, c("theta", "truth")], belt_truth[, c("theta", "truth")]), cfg_belt
    )),
    tar_target(ci_extra_obs, extra_cis(obs_result, sim_results, truth_table, belt)),
    tar_target(ci_methods, compare_ci_methods(sim_results, truth_table, belt)),
    tar_target(table_ci_methods, make_ci_table(ci_methods)),
    tar_target(fig_belt, save_belt_figure(belt, file.path(out_dir, "paper", ds)), format = "file"),
    tar_target(fig_ci_methods, save_ci_figure(ci_methods, file.path(out_dir, "paper", ds)),
               format = "file")
  )
} else {
  list()
}

# One target chain per dataset (static branching)
per_dataset <- tar_map(
  values = list(ds = active, hr_only = unname(is_hr)),
  names = "ds",

  # Data and observed analysis (paper Step 1)
  tar_target(base_file, cfg$datasets[[ds]]$file, format = "file"),
  tar_target(base_data, prepare_base_data(base_file, cfg$datasets[[ds]], cfg$datasets[[ds]]$base)),
  tar_target(obs_result, run_pipeline(base_data, cfg$datasets[[ds]]$pipeline)),

  # Simulated world: prediction model, inverse calibration of theta* (paper 2.3, 2.7.1), truths
  tar_target(dgm, fit_dgm(base_data, cfg$datasets[[ds]][["dgm"]], cfg$datasets[[ds]]$censoring)),
  tar_target(calibration, calibrate_observed(dgm, base_data, obs_result, crn_seeds,
                                             cfg$datasets[[ds]]$pipeline, cfg,
                                             cfg$datasets[[ds]]$calib_slope,
                                             cfg$datasets[[ds]]$tol_scale)),
  tar_target(theta_sim, c(calibration$theta_cal, 0)),
  tar_target(
    truth,
    compute_truth(dgm, base_data, theta_sim, cfg$datasets[[ds]]$pipeline, cfg),
    pattern = map(theta_sim)
  ),
  tar_target(truth_table, make_truth_table(truth, ds)),

  # Replicates: one branch per (scenario, batch); CRN seeds tied to replicate index
  tar_target(
    sim_results,
    run_batch(base_data, dgm, theta_sim, crn_seeds, batch_id, cfg$datasets[[ds]]$pipeline,
              cfg$batch_size)[, dataset := ds],
    pattern = cross(theta_sim, batch_id)
  ),
  # Independent replicates that check the adjusted analysis (own seed stream)
  tar_target(
    val_results,
    run_batch(base_data, dgm, theta_sim, val_seeds, val_batch_id, cfg$datasets[[ds]]$pipeline,
              cfg$batch_size)[, dataset := ds],
    pattern = cross(theta_sim, val_batch_id)
  ),

  # Performance and simple adjustment (paper 2.6-2.7)
  tar_target(performance, summarise_performance(
    sim_results, truth_table,
    favours_low = endpoint_info(cfg$datasets[[ds]]$pipeline)$favours_low
  )),
  tar_target(performance_rsimsum, rsimsum_check(sim_results, truth_table)),
  tar_target(adjustment, adjust_simple(
    obs_result, sim_results, truth_table, ds, cfg$datasets[[ds]]$pipeline$conf_level,
    endpoint_info(cfg$datasets[[ds]]$pipeline)$favours_low
  )),
  tar_target(validation, validate_adjustment(adjustment, val_results, truth_table,
                                             cfg$datasets[[ds]]$pipeline$conf_level)),

  # Nested simulation (paper 2.7.2): outer datasets at theta* and with no effect
  tar_target(
    nested_obs,
    run_nested_batch(base_data, dgm, calibration, val_seeds, nested_batch_id, cfg$datasets[[ds]],
                     cfg),
    pattern = map(nested_batch_id)
  ),
  tar_target(
    nested_null,
    run_nested_batch(base_data, dgm, NULL, val_seeds, nested_batch_id, cfg$datasets[[ds]], cfg),
    pattern = map(nested_batch_id)
  ),
  tar_target(nested, summarise_nested(nested_obs, nested_null, truth_table, adjustment,
                                      cfg$datasets[[ds]]$pipeline$conf_level,
                                      isTRUE(cfg$nested$shrink))),
  tar_target(perf_compare, compare_performance(
    sim_results, truth_table, validation, nested, cfg$datasets[[ds]]$pipeline$conf_level,
    endpoint_info(cfg$datasets[[ds]]$pipeline)$favours_low
  )),

  # Pilot study for the outer / inner split (paper Appendix): few outer, few inner runs, own seeds
  tar_target(
    pilot_nested,
    run_pilot_batch(base_data, dgm, calibration, pilot_seeds, pilot_batch_id, cfg$datasets[[ds]],
                    cfg),
    pattern = map(pilot_batch_id)
  ),
  tar_target(pilot_timing, time_inner_run(base_data, dgm, calibration$theta_cal,
                                          cfg$datasets[[ds]]$pipeline, cfg$pilot$timing_reps)),
  tar_target(allocation, recommend_allocation(pilot_nested, cfg$pilot, cfg$nested, pilot_timing)),
  tar_target(pilot_anova, anova_by_inner_runs(pilot_nested, c(5, 10, 20, 30, 50, 100, 250))),

  # Target for a pipeline that is not consistent (paper 4.6): unweighted hypothetical trial, same
  # replicates and inner runs, truth replaced (HR analysis only; other analyses give NULL)
  tar_target(truth_rct, if (hr_only) {
    truth_rct_grid(dgm, base_data, theta_sim, cfg$datasets[[ds]]$pipeline, cfg)
  }),
  tar_target(truth_table_rct, if (hr_only) make_truth_table(truth_rct, ds)),
  tar_target(adjustment_rct, if (hr_only) {
    adjust_simple(obs_result, sim_results, truth_table_rct, ds,
                  cfg$datasets[[ds]]$pipeline$conf_level, TRUE)
  }),
  tar_target(validation_rct, if (hr_only) {
    validate_adjustment(adjustment_rct, val_results, truth_table_rct,
                        cfg$datasets[[ds]]$pipeline$conf_level)
  }),
  tar_target(
    nested_obs_rct,
    if (hr_only) {
      rct_nested_update(nested_obs, base_data, dgm, calibration, val_seeds, cfg$datasets[[ds]], cfg,
                        cfg$datasets[[ds]]$pipeline$conf_level)
    },
    pattern = map(nested_obs)
  ),
  tar_target(nested_rct, if (hr_only) {
    summarise_nested(nested_obs_rct, nested_null, truth_table_rct, adjustment_rct,
                     cfg$datasets[[ds]]$pipeline$conf_level, isTRUE(cfg$nested$shrink))
  }),
  tar_target(perf_compare_rct, if (hr_only) {
    compare_performance(sim_results, truth_table_rct, validation_rct, nested_rct,
                        cfg$datasets[[ds]]$pipeline$conf_level, TRUE)
  }),
  tar_target(table_results_rct, if (hr_only) {
    make_results_table(obs_result, dgm, calibration, truth_table_rct, performance, adjustment_rct,
                       perf_compare_rct, nested_rct, cfg$datasets[[ds]])
  }),

  # Publication outputs (paper Section 4): tables and figures
  tar_target(raw_data, read_raw_data(base_file, cfg$datasets[[ds]]$base)),
  tar_target(table_baseline, make_baseline_table(raw_data, cfg$datasets[[ds]]$base,
                                                 cfg$datasets[[ds]])),
  tar_target(table_results, make_results_table(
    obs_result, dgm, calibration, truth_table, performance, adjustment, perf_compare, nested,
    cfg$datasets[[ds]]
  )),
  tar_target(paper_figures, save_paper_figures(
    perf_compare, sim_results, nested_obs, base_data, dgm, adjustment, cfg$datasets[[ds]],
    file.path(out_dir, "paper", ds)
  ), format = "file"),
  tar_target(balance, ccw_balance(base_data, cfg$datasets[[ds]]$pipeline,
                                   cfg$datasets[[ds]][["dgm"]]$covariates)),
  tar_target(fig_balance, save_balance_figure(balance, cfg$datasets[[ds]],
                                              file.path(out_dir, "paper", ds)), format = "file"),
  tar_target(fig_protocol, save_plot(
    plot_protocol(cfg$datasets[[ds]]$base, cfg$datasets[[ds]]$pipeline,
                  if (cfg$datasets[[ds]]$base == "maringe") 365.24 else cfg$datasets[[ds]]$horizon),
    file.path(out_dir, paste0("fig_protocol_", ds, ".png")), width = 8.5,
    height = if (cfg$datasets[[ds]]$base == "maringe") 4.6 else 5
  ), format = "file"),
  extra_targets
)

list(
  tar_target(crn_seeds, make_crn_seeds(cfg$seed_master, cfg$B)),
  tar_target(val_seeds, make_crn_seeds(cfg$seed_validation, cfg$B_validation)),
  tar_target(batch_id, seq_len(ceiling(cfg$B / cfg$batch_size))),
  tar_target(val_batch_id, seq_len(ceiling(cfg$B_validation / cfg$batch_size))),
  tar_target(nested_batch_id, seq_len(ceiling(cfg$nested$outer / cfg$nested$batch_size))),
  tar_target(pilot_seeds, make_crn_seeds(cfg$seed_pilot, cfg$pilot$outer)),
  tar_target(pilot_batch_id, seq_len(ceiling(cfg$pilot$outer / cfg$pilot$batch_size))),
  if (extras) tar_target(belt_batch_id, seq_len(ceiling(min(cfg_belt$B, cfg$B) / cfg$batch_size))),
  per_dataset,

  # Across datasets
  tar_combine(performance_all, per_dataset[["performance"]],
              command = data.table::rbindlist(list(!!!.x))),
  tar_combine(adjustment_all, per_dataset[["adjustment"]],
              command = data.table::rbindlist(list(!!!.x))),
  tar_combine(sim_all, per_dataset[["sim_results"]],
              command = data.table::rbindlist(list(!!!.x))),
  tar_combine(rsimsum_all, per_dataset[["performance_rsimsum"]],
              command = data.table::rbindlist(list(!!!.x))),
  tar_combine(calibration_all, per_dataset[["calibration"]], command = list(!!!.x)),
  tar_combine(dgm_all, per_dataset[["dgm"]], command = list(!!!.x)),
  tar_combine(obs_all, per_dataset[["obs_result"]], command = list(!!!.x)),
  tar_combine(base_all, per_dataset[["base_data"]], command = list(!!!.x)),

  # Outputs
  tar_target(fig_estimates, save_plot(plot_estimates(sim_all, performance_all, ds_labels),
                                      file.path(out_dir, "fig_estimates.png"), height = 5.5),
             format = "file"),
  tar_target(session, sessioninfo::session_info()),
  # How the run was parallelised (reported with the run time)
  tar_target(run_settings, record_run_settings(sim_all)),

  # Two-dataset comparison figure (two_datasets profile only)
  if (length(ds_labels) > 1) {
    tar_target(fig_comparison, save_plot(plot_comparison(performance_all, ds_labels),
                                         file.path(out_dir, "fig_comparison.png"), height = 5.5),
               format = "file")
  },
  # Colleague report on the Maringe case study (needs the Quarto CLI)
  if (do_report) {
    list(
      tar_target(report_qmd, "reports/qba_report.qmd", format = "file"),
      tar_target(report, render_report(
        report_qmd, targets::tar_path_store(), out_dir,
        deps = list(
          performance_maringe, adjustment_maringe, validation_maringe, nested_maringe,
          nested_obs_maringe, nested_null_maringe, perf_compare_maringe, calibration_maringe,
          dgm_maringe, obs_result_maringe, base_data_maringe, sim_results_maringe,
          truth_table_maringe, table_baseline_maringe, table_results_maringe,
          performance_rsimsum_maringe, balance_maringe, paper_figures_maringe,
          fig_protocol_maringe, fig_balance_maringe, run_settings, session
        )
      ), format = "file")
    )
  }
)
