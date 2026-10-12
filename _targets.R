# _targets.R - modelling-bias QBA simulation (CCW): observed-effect and no-effect scenarios
# for the datasets in cfg$active_datasets (default: Maringe only, as in the paper).
# Run with targets::tar_make(); smoke run: Sys.setenv(R_CONFIG_ACTIVE = "smoke");
# two-dataset comparison with report: Sys.setenv(R_CONFIG_ACTIVE = "two_datasets").
library(targets)
library(tarchetypes)

tar_source("R")
cfg <- config::get(file = "config/config.yml")
active <- unlist(cfg$active_datasets)
ds_labels <- vapply(cfg$datasets[active], `[[`, character(1), "label")
# Parallel workers need crew; without it the pipeline runs sequentially
use_crew <- cfg$workers > 1 && requireNamespace("crew", quietly = TRUE)

tar_option_set(
  packages = c("data.table", "survival", "checkmate"),
  seed = cfg$seed_master,
  controller = if (use_crew) crew::crew_controller_local(workers = cfg$workers) else NULL
)

# One target chain per dataset (static branching)
per_dataset <- tar_map(
  values = list(ds = active),
  names = "ds",

  # Data and observed analysis (paper Step 1)
  tar_target(base_file, cfg$datasets[[ds]]$file, format = "file"),
  tar_target(base_data, prepare_base_data(base_file, cfg$datasets[[ds]], ds)),
  tar_target(obs_result, run_pipeline(base_data, cfg$datasets[[ds]]$pipeline)),

  # Simulated world: reference DGM, calibration of theta* (paper Step 3), true values
  tar_target(dgm, fit_dgm(base_data, cfg$datasets[[ds]][["dgm"]], cfg$datasets[[ds]]$censoring)),
  tar_target(calibration, calibrate_theta(dgm, base_data, obs_result, cfg$datasets[[ds]]$pipeline, cfg)),
  tar_target(theta_sim, c(calibration$theta_star, 0)),
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

  # Summaries and correction (paper Steps 6-7)
  tar_target(performance, summarise_performance(sim_results, truth_table)),
  tar_target(performance_rsimsum, rsimsum_check(sim_results, truth_table)),
  tar_target(correction, correct_estimate(obs_result, performance, sim_results, ds)),

  # Publication outputs (paper Section 4): tables, adjusted performance, figures
  tar_target(raw_data, read_raw_data(base_file, ds)),
  tar_target(table_baseline, make_baseline_table(raw_data, ds, cfg$datasets[[ds]])),
  tar_target(perf_adjusted, evaluate_adjustment(sim_results, truth_table)),
  tar_target(table_results, make_results_table(
    obs_result, dgm, calibration, truth_table, performance, correction, perf_adjusted,
    cfg$datasets[[ds]]
  )),
  tar_target(paper_figures, save_paper_figures(
    performance, perf_adjusted, sim_results, base_data, dgm, correction$z_crit_calibrated,
    cfg$datasets[[ds]], file.path("reports", "paper", ds)
  ), format = "file"),
  tar_target(fig_protocol, save_plot(
    plot_protocol(ds, cfg$datasets[[ds]]$pipeline,
                  if (ds == "maringe") 365.24 else cfg$datasets[[ds]]$horizon),
    file.path("reports", paste0("fig_protocol_", ds, ".png")), width = 8.5,
    height = if (ds == "maringe") 4.6 else 5
  ), format = "file")
)

list(
  tar_target(crn_seeds, make_crn_seeds(cfg$seed_master, cfg$B)),
  tar_target(batch_id, seq_len(ceiling(cfg$B / cfg$batch_size))),
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
                                      "reports/fig_estimates.png", height = 5.5), format = "file"),
  tar_target(session, sessioninfo::session_info()),

  # Two-dataset comparison and colleague report (two_datasets profile only)
  if (length(active) > 1) {
    tar_target(fig_comparison, save_plot(plot_comparison(performance_all, ds_labels),
                                         "reports/fig_comparison.png", height = 5.5),
               format = "file")
  },
  if (isTRUE(cfg$render_report)) tar_quarto(report, "reports/qba_report.qmd", quiet = FALSE)
)
