# _targets.R - modelling-bias QBA simulation (CCW): observed-effect and no-effect scenarios
# for the datasets in cfg$active_datasets (default: Maringe only, as in the paper).
# Run with targets::tar_make(); smoke run: Sys.setenv(R_CONFIG_ACTIVE = "smoke");
# two-dataset comparison: Sys.setenv(R_CONFIG_ACTIVE = "two_datasets").
library(targets)
library(tarchetypes)

tar_source("R")
cfg <- config::get(file = "config/config.yml")
# config attaches the absolute path of the file; dropping it keeps the store valid when the
# project is moved or cloned elsewhere (every target depends on cfg)
attributes(cfg) <- list(names = names(cfg))
# Confidence-belt settings live in their own file so that adding them did not invalidate cfg
cfg_belt <- config::get(file = "config/belt.yml")
attributes(cfg_belt) <- list(names = names(cfg_belt))
active <- unlist(cfg$active_datasets)
ds_labels <- vapply(cfg$datasets[active], `[[`, character(1), "label")
# Outputs of the default (paper) profile go to reports/; other profiles get their own folder so a
# precision run cannot overwrite the paper figures
profile <- Sys.getenv("R_CONFIG_ACTIVE", "default")
out_dir <- if (profile == "default") "reports" else file.path("reports", profile)
# Report on the Maringe case study (into out_dir): not in smoke runs, and only if Quarto is installed.
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

# One target chain per dataset (static branching)
per_dataset <- tar_map(
  values = list(ds = active),
  names = "ds",

  # Data and observed analysis (paper Step 1)
  tar_target(base_file, cfg$datasets[[ds]]$file, format = "file"),
  tar_target(base_data, prepare_base_data(base_file, cfg$datasets[[ds]], ds)),
  tar_target(obs_result, run_pipeline(base_data, cfg$datasets[[ds]]$pipeline)),

  # Simulated world: reference DGM, calibration of theta* (paper Step 3), true values.
  # One branch per large dataset, so calibration steps and truths run in parallel; every
  # (theta, seed) evaluation is reused wherever it is needed again (CRN)
  tar_target(dgm, fit_dgm(base_data, cfg$datasets[[ds]][["dgm"]], cfg$datasets[[ds]]$censoring)),
  tar_target(calib_theta_1, calib_start(dgm)),
  tar_target(
    calib_rows_1,
    truth_point(dgm, base_data, calib_theta_1, calib_seed_id, truth_seeds,
                cfg$datasets[[ds]]$pipeline, cfg),
    pattern = cross(calib_theta_1, calib_seed_id)
  ),
  # cfg[["calibration"]], not cfg$calibration: tar_map() renames the symbol `calibration` to this
  # dataset's calibration target (cf. [["dgm"]] above), which would pass NULL
  tar_target(calib_theta_2, calib_next(calib_rows_1, obs_result$estimate, cfg[["calibration"]])),
  tar_target(
    calib_rows_2,
    truth_point(dgm, base_data, calib_theta_2, calib_seed_id, truth_seeds,
                cfg$datasets[[ds]]$pipeline, cfg, known = calib_rows_1),
    pattern = map(calib_seed_id)
  ),
  tar_target(calib_theta_3, calib_next(rbind(calib_rows_1, calib_rows_2), obs_result$estimate,
                                       cfg[["calibration"]])),
  tar_target(
    calib_rows_3,
    truth_point(dgm, base_data, calib_theta_3, calib_seed_id, truth_seeds,
                cfg$datasets[[ds]]$pipeline, cfg, known = rbind(calib_rows_1, calib_rows_2)),
    pattern = map(calib_seed_id)
  ),
  # Serial secant steps only if the two parallel stages did not reach the tolerance
  tar_target(calibration, finish_calibration(
    rbind(calib_rows_1, calib_rows_2, calib_rows_3), dgm, base_data, obs_result,
    cfg$datasets[[ds]]$pipeline, cfg, truth_seeds
  )),
  tar_target(theta_sim, c(calibration$theta_star, 0)),
  # Calibration already evaluated theta* and 0 on the first calibration seeds
  tar_target(
    truth_rows,
    truth_point(dgm, base_data, theta_sim, truth_seed_id, truth_seeds,
                cfg$datasets[[ds]]$pipeline, cfg, known = calibration$evals),
    pattern = cross(theta_sim, truth_seed_id)
  ),
  tar_target(truth, summarise_truth(truth_rows)),
  tar_target(truth_table, make_truth_table(truth, ds)),

  # Replicates: one branch per (scenario, batch); CRN seeds tied to replicate index
  tar_target(
    sim_results,
    run_batch(base_data, dgm, theta_sim, crn_seeds, batch_id, cfg$datasets[[ds]]$pipeline,
              cfg$batch_size)[, dataset := ds],
    pattern = cross(theta_sim, batch_id)
  ),

  # Summaries and correction (paper Steps 6-7)
  tar_target(performance, summarise_performance(sim_results, truth_table)),
  tar_target(performance_rsimsum, rsimsum_check(sim_results, truth_table)),
  # Confidence belt: Chebyshev-Lobatto nodes over the psi range where CI endpoints fall, mapped
  # to theta through the psi(theta) points already evaluated; truths and replicates (same CRN
  # seeds) at each node; replicates topped up until the quantiles of T reach the target MCSE
  tar_target(belt_range, belt_psi_range(obs_result, sim_results, cfg_belt$tail, cfg_belt$expand)),
  tar_target(belt_theta, belt_theta_plan(
    belt_range, truth_table, rbind(calibration$evals, truth_rows), cfg_belt
  )),
  tar_target(
    belt_truth_rows,
    truth_point(dgm, base_data, belt_theta, belt_seed_id, truth_seeds,
                cfg$datasets[[ds]]$pipeline, cfg, known = calibration$evals),
    pattern = cross(belt_theta, belt_seed_id)
  ),
  tar_target(belt_truth, summarise_truth(belt_truth_rows)),
  tar_target(
    belt_sims_1,
    run_batch(base_data, dgm, belt_theta, belt_seeds, belt_batch_id, cfg$datasets[[ds]]$pipeline,
              cfg$batch_size)[, dataset := ds],
    pattern = cross(belt_theta, belt_batch_id)
  ),
  tar_target(belt_plan, belt_topup_plan(belt_sims_1, belt_truth, cfg_belt, cfg$batch_size)),
  tar_target(
    belt_sims_2,
    run_belt_topup(base_data, dgm, belt_theta, belt_seeds, belt_topup_id, belt_plan,
                   cfg$datasets[[ds]]$pipeline, cfg$batch_size),
    pattern = cross(belt_theta, belt_topup_id)
  ),
  tar_target(belt, fit_belt(
    rbind(sim_results, belt_sims_1, belt_sims_2, fill = TRUE),
    rbind(truth_table[, .(theta, truth)], belt_truth[, .(theta, truth)]), cfg_belt
  )),
  tar_target(correction, correct_estimate(obs_result, performance, sim_results, ds, belt = belt)),

  # Publication outputs (paper Section 4): tables, adjusted performance, figures
  tar_target(raw_data, read_raw_data(base_file, ds)),
  tar_target(table_baseline, make_baseline_table(raw_data, ds, cfg$datasets[[ds]])),
  tar_target(perf_adjusted, evaluate_adjustment(sim_results, truth_table)),
  # Alternative CIs (Wald, QBA widening, bootstrap-t, confidence belt): coverage, table, figures
  tar_target(ci_methods, compare_ci_methods(sim_results, truth_table, belt)),
  tar_target(fig_belt, save_belt_figure(belt, file.path(out_dir, "paper", ds)), format = "file"),
  tar_target(table_ci_methods, make_ci_table(ci_methods)),
  tar_target(fig_ci_methods, save_ci_figure(ci_methods, file.path(out_dir, "paper", ds)),
             format = "file"),
  tar_target(table_results, make_results_table(
    obs_result, dgm, calibration, truth_table, performance, correction, perf_adjusted,
    cfg$datasets[[ds]]
  )),
  tar_target(paper_figures, save_paper_figures(
    performance, perf_adjusted, sim_results, base_data, dgm, correction$z_crit_calibrated,
    cfg$datasets[[ds]], file.path(out_dir, "paper", ds)
  ), format = "file"),
  tar_target(balance, ccw_balance(base_data, cfg$datasets[[ds]]$pipeline,
                                   cfg$datasets[[ds]][["dgm"]]$covariates)),
  tar_target(fig_balance, save_balance_figure(balance, cfg$datasets[[ds]],
                                              file.path(out_dir, "paper", ds)), format = "file"),
  tar_target(fig_protocol, save_plot(
    plot_protocol(ds, cfg$datasets[[ds]]$pipeline,
                  if (ds == "maringe") 365.24 else cfg$datasets[[ds]]$horizon),
    file.path(out_dir, paste0("fig_protocol_", ds, ".png")), width = 8.5,
    height = if (ds == "maringe") 4.6 else 5
  ), format = "file")
)

list(
  tar_target(crn_seeds, make_crn_seeds(cfg$seed_master, cfg$B)),
  tar_target(batch_id, seq_len(ceiling(cfg$B / cfg$batch_size))),
  # Truth seeds: CRN across theta; calibration, main truths and belt truths use prefixes of them
  tar_target(truth_seeds, make_crn_seeds(
    cfg$seed_truth, max(cfg$truth$n_rep, cfg$calibration$n_rep, cfg_belt$truth_n_rep)
  )),
  tar_target(calib_seed_id, seq_len(cfg$calibration$n_rep)),
  tar_target(truth_seed_id, seq_len(cfg$truth$n_rep)),
  tar_target(belt_seed_id, seq_len(cfg_belt$truth_n_rep)),
  # Belt replicates: first stage B_start per node, top-up batches up to B_max where needed
  tar_target(belt_seeds, crn_seeds[seq_len(min(cfg_belt$B_max, cfg$B))]),
  tar_target(belt_batch_id, seq_len(ceiling(min(cfg_belt$B_start, cfg$B) / cfg$batch_size))),
  tar_target(belt_topup_id, {
    n1 <- length(belt_batch_id)
    # At least one id so the pattern is never empty; unneeded batches return no rows
    seq.int(n1 + 1L, max(n1 + 1L, ceiling(min(cfg_belt$B_max, cfg$B) / cfg$batch_size)))
  }),
  per_dataset,

  # Across datasets
  tar_combine(performance_all, per_dataset[["performance"]],
              command = data.table::rbindlist(list(!!!.x))),
  tar_combine(correction_all, per_dataset[["correction"]],
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
  if (length(active) > 1) {
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
          performance_maringe, correction_maringe, calibration_maringe, dgm_maringe,
          belt_plan_maringe,
          obs_result_maringe, base_data_maringe, sim_results_maringe, truth_table_maringe,
          perf_adjusted_maringe, table_baseline_maringe, table_results_maringe,
          ci_methods_maringe, table_ci_methods_maringe, fig_ci_methods_maringe, belt_maringe,
          fig_belt_maringe,
          performance_rsimsum_maringe, balance_maringe, paper_figures_maringe,
          fig_protocol_maringe, fig_balance_maringe, run_settings, session
        )
      ), format = "file")
    )
  }
)
