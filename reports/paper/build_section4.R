# Word draft of paper Section 4 (Methods text, Results text, tables, figures) from the
# targets store. Run from the project root after targets::tar_make():
#   source("reports/paper/build_section4.R")
# Kept outside the pipeline while the wording is being iterated; reads run times from tar_meta().
# Numbering follows manuscript v0.65 (24 Oct 2026): Table 1, Figures 2-3 in 4.1-4.4; 4.6 and the
# Appendix are added from the extra analyses when their targets exist in the store.

# Text contains en dashes and Greek letters; R >= 4.2 on Windows is UTF-8 by default
if (!isTRUE(l10n_info()[["UTF-8"]])) message("Run R in a UTF-8 locale, e.g. LANG=C.UTF-8")
library(data.table)
library(officer)
library(flextable)
cfg <- config::get(file = "config/config.yml")
# Same profile as the run: results from _targets (paper run) or _targets_<profile>, written
# to reports/ or reports/<profile>/
profile <- Sys.getenv("R_CONFIG_ACTIVE", "default")
store <- if (profile == "default") "_targets" else paste0("_targets_", profile)
out_root <- if (profile == "default") "reports" else file.path("reports", profile)
ds <- "maringe"
stopifnot(ds %in% unlist(cfg$active_datasets))
# Targets of the paper dataset under short names (e.g. performance_maringe -> performance)
for (nm in c("table_baseline", "table_results", "perf_compare", "performance", "adjustment",
             "nested", "calibration", "dgm", "obs_result", "base_data", "truth_table",
             "sim_results", "balance")) {
  assign(nm, targets::tar_read_raw(paste0(nm, "_", ds), store = store))
}
spec <- cfg$datasets[[ds]]
fig_dir <- file.path(out_root, "paper", ds)
out_file <- file.path(out_root, "paper", "section4_draft.docx")

# ---- Numbers used in the text --------------------------------------------------
f1 <- function(x) formatC(x, format = "f", digits = 1)
f2 <- function(x) formatC(x, format = "f", digits = 2)
f3 <- function(x) formatC(x, format = "f", digits = 3)
f4 <- function(x) formatC(x, format = "f", digits = 4)
pc <- function(x, d = 1) paste0(formatC(100 * x, format = "f", digits = d), "%")
n0 <- function(x) format(x, big.mark = ",")
# HRs and CIs with 4 decimals; type I error with 3 (percent)
hr_ci <- function(e, lo, hi) paste0(f4(exp(e)), " (95% CI ", f4(exp(lo)), "–", f4(exp(hi)), ")")
po <- performance[scenario == "observed"]
pn <- performance[scenario == "null"]
adj <- adjustment
ns <- nested$summary
cal <- calibration
bd <- base_data
tn <- truth_table[scenario == "null"]
P <- function(an, meth, m) perf_compare[analysis == an & method == meth & measure == m]
sa <- function(m) P("Simple adjustment", "Non-parametric", m)
na <- function(m) P("Nested adjustment", "Non-parametric", m)
mc_ci <- function(v, mc) paste0(f3(v - 1.96 * mc), " to ", f3(v + 1.96 * mc))
off <- function(r) abs(r$value - r$ref) > 1.96 * r$mcse
lim_sd <- function(n) sqrt(2 * 0.025 * 0.975 / n)
z_tot <- function(r, n) (r$value - r$ref) / sqrt(r$mcse^2 + lim_sd(n)^2)
B <- po$B
n_val <- sa("95% CI coverage")$n
n_fail <- po$n_fail + pn$n_fail
n_trt <- sum(is.finite(bd$trt_time))
# Patterns repeat the time of their branches, so only stems and branches are summed
source("R/utils.R", local = TRUE)
# Only the current pipeline: the store also keeps branch records of earlier runs
meta <- targets::tar_meta(fields = c(name, type, time, seconds, children), targets_only = TRUE,
                          store = store)
meta <- current_meta(meta)
secs_of <- function(prefix) sum(meta[grepl(paste0("^", prefix, "_", ds, "_"), meta$name), ]$seconds, na.rm = TRUE)
secs_sim <- secs_of("sim_results")
secs_nested <- secs_of("nested_obs") + secs_of("nested_null")
secs_all <- sum(meta$seconds, na.rm = TRUE)
heavy <- meta[grepl(paste0("^(calibration|truth_rows|sim_results|val_results|nested_obs|nested_null)_", ds), meta$name), ]
mins_wall <- as.numeric(difftime(max(heavy$time), min(heavy$time - heavy$seconds), units = "mins"))
rs <- targets::tar_read(run_settings, store = store)
# Workers actually used, counted from the timestamps of the simulation branches
rs$workers <- observed_workers(meta)
# R version of the run itself (session target), not of the session building this draft
r_ver <- sub("^R version ([0-9.]+).*", "\\1", targets::tar_read(session, store = store)$platform$version)
n_surg_day1 <- sum(bd$trt_time <= 1)
n_big <- function(n_target) n0(ceiling(n_target / nrow(bd)) * nrow(bd))

# ---- Document helpers ----------------------------------------------------------
fp_body <- fp_text(font.family = "Calibri", font.size = 11)
fp_bold <- update(fp_body, bold = TRUE)
fp_note <- fp_text(font.family = "Calibri", font.size = 9.5, italic = TRUE, color = "#555555")
para <- function(doc, ...) body_add_fpar(doc, fpar(..., fp_p = fp_par(padding.bottom = 6)))
txt <- function(doc, x) para(doc, ftext(x, fp_body))
note <- function(doc, x) para(doc, ftext(x, fp_note))
# Plain bold headings: the default Word template numbers its heading styles
fp_h1 <- update(fp_body, bold = TRUE, font.size = 14)
fp_h2 <- update(fp_body, bold = TRUE, font.size = 12)
head1 <- function(doc, x) body_add_fpar(doc, fpar(ftext(x, fp_h1), fp_p = fp_par(padding.bottom = 6)))
head2 <- function(doc, x) {
  body_add_fpar(doc, fpar(ftext(x, fp_h2), fp_p = fp_par(padding.top = 10, padding.bottom = 4)))
}
caption <- function(doc, label, x) para(doc, ftext(label, fp_bold), ftext(x, fp_body))

# Tables in the layout of Table S2: no vertical lines, rules under the header and at the end
style_table <- function(ft) {
  ft |>
    font(fontname = "Calibri", part = "all") |>
    fontsize(size = 10, part = "all") |>
    border_remove() |>
    hline_top(part = "header", border = fp_border(width = 1)) |>
    hline_bottom(part = "header", border = fp_border(width = 0.75)) |>
    hline_bottom(part = "body", border = fp_border(width = 1)) |>
    padding(padding.top = 1, padding.bottom = 1, part = "all") |>
    line_spacing(space = 1, part = "all") |>
    bold(part = "header")
}

baseline_ft <- function(tb) {
  d <- tb[, .(label, n1, p1, n0, p0, smd)]
  ft <- flextable(d) |>
    set_header_labels(label = "", n1 = "N", p1 = "%", n0 = "N", p0 = "%", smd = "") |>
    add_header_row(values = c("", attr(tb, "groups"), "Standardised mean differences (%)"),
                   colwidths = c(1, 2, 2, 1)) |>
    merge_at(i = 1:2, j = 6, part = "header") |>
    align(j = 2:6, align = "center", part = "all") |>
    style_table()
  # Section heads and Total in bold, levels indented as in Table S2
  bold_rows <- which(tb$indent == 0L)
  ft <- bold(ft, i = bold_rows, j = 1, part = "body")
  ft <- padding(ft, i = which(tb$indent == 1L), j = 1, padding.left = 18, part = "body")
  # The "Median / Q1 Q3" sub-header row
  sub <- which(tb$n1 == "Median")
  ft <- bold(ft, i = sub, part = "body")
  width(ft, j = 1:6, width = c(2.4, 0.7, 0.8, 0.7, 0.8, 1.3))
}

two_col_ft <- function(d, h1, h2) {
  ft <- flextable(d) |> set_header_labels(values = setNames(c(h1, h2), names(d))) |> style_table()
  width(ft, j = 1:2, width = c(4.4, 2.3))
}

perf_ft <- function(d) {
  flextable(d) |>
    set_header_labels(measure = "Performance measure", ccw = "CCW as published",
                      simple = "Simple adjustment", nested = "Nested simulation") |>
    align(j = 2:4, align = "center", part = "all") |>
    style_table() |>
    width(j = 1:4, width = c(2.9, 1.35, 1.35, 1.5))
}

ci_var_ft <- function(d) {
  flextable(d) |>
    set_header_labels(analysis = "Analysis", method = "CI", coverage = "Coverage, % (MCSE)",
                      below = "Below truth, %", above = "Above truth, %") |>
    merge_v(j = "analysis") |>
    valign(j = 1, valign = "top") |>
    align(j = 3:5, align = "center", part = "all") |>
    style_table() |>
    width(j = 1:5, width = c(1.5, 2.2, 1.3, 0.9, 0.9))
}

dgm_ft <- function(d) {
  flextable(d[, .(model, term, hr, ci)]) |>
    set_header_labels(model = "Model", term = "Term", hr = "HR", ci = "95% CI") |>
    merge_v(j = "model") |>
    valign(j = 1, valign = "top") |>
    style_table() |>
    width(j = 1:4, width = c(1.8, 2.3, 0.8, 1.4))
}

add_fig <- function(doc, file, w, h) body_add_img(doc, src = file, width = w, height = h)

# ---- Build ---------------------------------------------------------------------
doc <- read_docx()
doc <- head1(doc, "Section 4 – Example (draft for insertion)")
doc <- note(doc, paste0(
  "Generated by reports/paper/build_section4.R from the results store ", store, " (run of ",
  format(Sys.Date(), "%d %b %Y"), "; ", n0(B), " simulated studies per scenario, ", n0(n_val),
  " independent studies per scenario for the check of the adjustment, nested simulation with ",
  ns$n_outer, " outer datasets × ", ns$n_inner, " inner studies; protocol v0.4). ",
  "Text in square brackets needs a decision or a check by the authors. Table and figure numbers ",
  "follow manuscript v0.64."
))

# 4.1 Used dataset
doc <- head2(doc, "4.1 Used dataset")
doc <- txt(doc, paste0(
  "We reanalyse the example dataset of Maringe et al. [35] to check whether the applied clone-censor-",
  "weighting (CCW) approach is unbiased, has nominal 95% CI coverage and no type I error inflation. ",
  "The dataset comprises ", n0(nrow(bd)), " simulated elderly patients with stage I or II lung cancer, ",
  "generated by the authors to mimic their cancer registry study and distributed with the tutorial ",
  "[it is the simulated supplementary dataset, not the registry data of their Table S2; data from the ccwr ",
  "R package]. Time zero is the date of diagnosis and patients are followed for one year. The ",
  "treatment of interest is surgery within 6 months of diagnosis, received by ", n0(n_trt), " patients (",
  pc(n_trt / nrow(bd)), "), ", n_surg_day1, " of them on day 1; the median time to surgery among operated ",
  "patients was ", stats::median(bd$trt_time[is.finite(bd$trt_time)]), " days. The endpoint is death ",
  "within one year, observed in ", n0(sum(bd$status)), " patients (", pc(mean(bd$status)), "). All ",
  "patients alive at one year were followed until day 365, so the only censoring is administrative ",
  "(Table 1)."
))
doc <- caption(doc, "Table 1: ", "Baseline characteristics")
doc <- body_add_flextable(doc, baseline_ft(table_baseline), align = "left")
doc <- note(doc, "Standardised mean differences: one value for binary characteristics, one per category for performance status and Charlson index, continuous scale for age and deprivation score. Deprivation score as coded in the simulated dataset (0–1 scale).")
doc <- caption(doc, "Figure 2: ", "Clone-censoring: how example patients are cloned into both strategies, when each clone is artificially censored, and which part of follow-up carries a censoring weight. Shaded area: grace period.")
doc <- add_fig(doc, file.path(out_root, paste0("fig_protocol_", ds, ".png")), 6.3, 6.3 * 4.6 / 8.5)

# 4.2 Model
doc <- head2(doc, "4.2 Model")
doc <- txt(doc, paste0(
  "Reference model. The reference model is the CCW analysis of Maringe et al. [35], comparing the ",
  "strategies 'surgery within 6 months of diagnosis' and 'no surgery within 6 months' (grace period ",
  spec$pipeline$grace, " days). Each patient is cloned into both strategies. Clones in the no-surgery ",
  "strategy are censored on the day of surgery if surgery happens within the grace period. Clones in ",
  "the surgery strategy are censored at the end of the grace period if they are alive and have not been ",
  "operated; patients who die within the grace period without surgery remain in both strategies. In the ",
  "no-surgery strategy, the probability of remaining uncensored is estimated from a Cox model for time to ",
  "surgery and the weights are updated at every surgery date; in the surgery strategy, a logistic model for ",
  "surgery within the grace period among patients alive at its end gives the weights from that time ",
  "onwards. Both models include age, sex, performance status, stage, deprivation score, Charlson index and ",
  "emergency presentation. The strategies are compared by a weighted Cox model with the strategy as the ",
  "only covariate, Efron ties and a robust, patient-clustered variance (Figure 2)."
))
doc <- txt(doc, paste0(
  "Preserved data. Not only the baseline data are maintained but also each patient's observed date of ",
  "surgery; patients who were not operated keep no surgery date. Surgery is counted as received only if ",
  "its date precedes the simulated end of follow-up."
))
doc <- txt(doc, paste0(
  "Prediction model. Times to death are predicted from a Weibull proportional hazards model fitted to the ",
  "original data by maximum likelihood, with the same seven covariates and surgery as a time-varying ",
  "covariate, h(t | X) = λρt^(ρ−1) exp(γ′X + θ·1{t ≥ s}), where s is the patient's surgery date, so the ",
  "effect of surgery acts only after surgery has taken place (Table S1). Death times are drawn by ",
  "inverting each patient's cumulative hazard; patients alive at day 365 are censored there."
))
doc <- txt(doc, paste0(
  "Simulation. ", n0(B), " datasets were simulated at the calibrated generating value θ* (Section 4.3) ",
  "and ", n0(B), " with no treatment effect (θ = 0), each reanalysed with the reference model. Replicate b ",
  "used the same random numbers in both scenarios (common random numbers). Performance measures were the ",
  "bias of the log HR, the empirical and the model-based (robust) standard error, the coverage of the ",
  "nominal 95% CI, the root mean squared error (RMSE) with its square (MSE) and the rejection rate of ",
  "H0: HR = 1 under no effect (type I error), each with its ",
  "Monte Carlo standard error (MCSE); the MCSE of the bias includes the Monte Carlo error of the true value. ",
  "Failed model fits were recorded and excluded (n = ", n_fail, ")."
))
doc <- txt(doc, paste0(
  "Adjustment (Section 2.7). All adjustments were non-parametric. The bias-adjusted estimate is ",
  "ψ(θ*) = ψ_obs − b(θ*). The adjusted 95% CI adds the 2.5th and 97.5th percentiles of the simulated ",
  "estimates at θ*, centred at their mean, to the bias-adjusted estimate. The adjusted test rejects H0 ",
  "when the observed estimate, minus the mean of the no-effect estimates, falls outside the 2.5th and ",
  "97.5th percentiles of the centred no-effect estimates; the simulation-calibrated p-value is twice the ",
  "smaller tail proportion. The adjusted analysis was checked in ", n0(n_val), " further, independent ",
  "datasets per scenario. In addition, a nested simulation (Section 2.7.2) was performed with ",
  ns$n_outer, " outer datasets: in each, the Weibull prediction model was refitted, θ was recalibrated ",
  "with ", ns$n_inner, " inner datasets per search step, the bias was re-estimated against the true value ",
  "under the refitted model, and the non-parametric CI was constructed from the inner datasets. Coverage ",
  "was the proportion of outer CIs containing ψ(θ*). With no effect, each of ", ns$n_null, " outer ",
  "datasets received its own rejection limits from ", ns$n_inner, " inner datasets under the refitted model."
))
doc <- txt(doc, paste0(
  "Software and run time. The simulation was programmed in R ", r_ver, " (packages survival, data.table ",
  "and targets; parallel workers via crew) and is available in the Online Supplement. It took ",
  f1(mins_wall), " minutes with ", rs$workers, " parallel worker", if (isTRUE(rs$workers > 1)) "s" else "",
  " (", f1(secs_all / 60), " minutes of computing time, of which ", f1(secs_sim / 60), " minutes for the ",
  "2 × ", n0(B), " simulated datasets, ", f2(secs_sim / (2 * B)), " seconds each, and ",
  f1(secs_nested / 60), " minutes for the nested simulation)."
))

# 4.3 Calibration
doc <- head2(doc, "4.3 Calibration")
doc <- txt(doc, paste0(
  "Calibration was performed as described in Section 2.7.1. The Weibull treatment parameter θ represents ",
  "the conditional effect of surgery once received, whereas the CCW analysis estimates a marginal contrast ",
  "between treatment strategies over the full follow-up period. These quantities differ because patients ",
  "assigned to the surgery strategy remain untreated until surgery, the strategy hazards are ",
  "non-proportional, and hazard ratios are non-collapsible. A secant search with ", n0(cal$B),
  " simulated datasets per step (common random numbers) identified the Weibull treatment parameter for ",
  "which the mean CCW estimate equalled the observed CCW estimate: a Weibull generating HR of ",
  f4(exp(cal$theta_cal)), " (θ* = ", f3(cal$theta_cal), "; fitted θ̂ = ", f3(dgm$theta_hat), ") reproduces ",
  "the observed CCW HR of ", f4(exp(obs_result$estimate)), " (mean over all ", n0(B), " datasets ",
  f4(exp(adj$mean_est)), "). The true CCW HR at θ*, defined as the probability limit of the reference ",
  "analysis and approximated by ", cfg$truth$n_rep, " simulated populations of ", n_big(cfg$truth$n_target),
  " patients, is ", f4(exp(adj$truth)), " (MCSE of the log HR ", f3(adj$truth_mcse), "). No calibration ",
  "was required for the type I error analysis because both the conditional and marginal treatment ",
  "effects equal HR = 1 under the null hypothesis. With unlimited data the CCW analysis would ",
  "nevertheless estimate HR ", f4(exp(tn$psi_large_sample)), " when surgery has no effect."
))

# 4.4 Goodness of fit
doc <- head2(doc, "4.4 Goodness of fit")
doc <- caption(doc, "Figure 3: ", "Fit of the Weibull prediction model: Kaplan–Meier estimate of the original data and the survival predicted by the model (averaged over patients, with each patient's observed surgery date).")
doc <- add_fig(doc, file.path(fig_dir, "fig_model_fit.png"), 4.6, 4.6 * 4.6 / 6)
doc <- txt(doc, "Because there was no loss to follow-up, there was no need to model the censoring distribution.")

# 4.5 Results
doc <- head2(doc, "4.5 Results")
bias_sig <- abs(po$bias) > 1.96 * po$bias_mcse_total
cov_sig <- abs(po$coverage - 0.95) > 1.96 * po$coverage_mcse
t1_sig <- abs(pn$reject_0 - 0.05) > 1.96 * pn$reject_0_mcse
doc <- txt(doc, paste0(
  "The Maringe example showed for the selected analysis pipeline ",
  if (bias_sig) "a statistically significant bias of " else "a bias of ", f3(po$bias),
  " on the log HR scale (95% MC interval ", mc_ci(po$bias, po$bias_mcse_total), "; ",
  f1(abs(100 * po$std_bias)), "% of the empirical SE), 95% CI coverage of ", pc(po$coverage),
  " (MCSE ", pc(po$coverage_mcse), ")", if (cov_sig) ", significantly below 95%," else "",
  " and a type I error rate of ", pc(pn$reject_0, 3), " (MCSE ", pc(pn$reject_0_mcse, 3), ")",
  if (t1_sig) ", significantly above 5%" else "", "; ", pc(pn$reject_0_fav, 3),
  " of the no-effect datasets were rejections in favour of surgery. The empirical SE was ", f3(po$emp_se),
  ", while the average robust SE was ", f3(po$mod_se), ". The RMSE was ", f3(po$rmse), " (MSE ",
  f4(po$mse), ", MCSE ", f4(po$mse_mcse), ") (Table 3)."
))
doc <- txt(doc, paste0(
  "The original value of HR = ", f4(exp(obs_result$estimate)), " was corrected to HR = ", f4(exp(adj$est_adj)),
  " (bias MCSE ", f3(adj$bias_mcse), " on the log scale). The non-parametric 95% CI of [",
  f4(exp(adj$ci_lo_np)), "; ", f4(exp(adj$ci_hi_np)), "] replaces the original [", f4(exp(obs_result$ci_lo)),
  "; ", f4(exp(obs_result$ci_hi)), "]. The adjusted 5% test rejects H0 if the CCW HR is below ",
  f4(exp(adj$null_mean + adj$null_q_lo)), " or above ", f4(exp(adj$null_mean + adj$null_q_hi)),
  "; the simulation-calibrated p-value is ", f3(adj$p_np), " (original p = ", f3(adj$p_obs), "). ",
  "In ", n0(n_val), " independent datasets per scenario, the adjusted analysis had a bias of ",
  f3(sa("Bias")$value), " (MCSE ", f3(sa("Bias")$mcse), "), 95% CI coverage of ",
  pc(sa("95% CI coverage")$value), " (MCSE ", pc(sa("95% CI coverage")$mcse), ") and a type I error ",
  "rate of ", pc(sa("Type I error")$value, 3), " (MCSE ", pc(sa("Type I error")$mcse, 3), ") and an RMSE of ",
  f3(sa("RMSE")$value), " (MSE ", f4(sa("MSE")$value), ")",
  # A learned rule's two 2.5% limits come from B studies, so its true rate varies by lim_sd
  if (off(sa("Type I error")) && abs(z_tot(sa("Type I error"), B)) <= 1.96) {
    paste0(", a deviation from 5% within the Monte Carlo error of the rejection limits themselves, ",
           "which were estimated from ", n0(B), " no-effect datasets (SD of the implied rejection ",
           "rate ", pc(lim_sd(B), 2), "; z = ", f2(z_tot(sa("Type I error"), B)), ")")
  } else "",
  " (Table 3, Figure 4)."
))
doc <- txt(doc, paste0(
  "Nested simulation. The calibrated generating value varied across the ", ns$n_outer_ok, " outer datasets ",
  "(SD of θ ", f3(ns$theta_r_sd), "); the estimated bias correction had mean ", f3(ns$bias_r_mean),
  " and SD ", f3(ns$bias_r_sd), " (inner Monte Carlo error ", f3(ns$bias_r_within_mcse), "). The outer bias-adjusted estimates had 95% CI coverage of ", pc(na("95% CI coverage")$value),
  " (MCSE ", pc(na("95% CI coverage")$mcse), "), and the outer rejection limits gave a type I error rate of ",
  pc(na("Type I error")$value, 3), " (MCSE ", pc(na("Type I error")$mcse, 3), "). Incorporating the ",
  "variation of the bias correction, the 95% CI of the bias-adjusted estimate is [", f4(exp(ns$ci_lo_np)),
  "; ", f4(exp(ns$ci_hi_np)), "] (Table 2, Figure 6)."
))
doc <- note(doc, paste0(
  "[Draft text to update: bias −0.015, type I error 7.5%/7.6%, coverage 93.6%, adjusted HR 0.6301 (0.37; 1.08), ",
  "critical value 2.17, Weibull HR 0.49. This run: bias ", f3(po$bias), ", type I error ", pc(pn$reject_0, 3),
  ", coverage ", pc(po$coverage), ", adjusted HR ", f4(exp(adj$est_adj)), ", Weibull HR ", f4(exp(cal$theta_cal)),
  ". The critical value is replaced by the non-parametric rejection limits.]"
))

doc <- caption(doc, "Table 2: ", "Results (1): observed analysis, calibration and QBA adjustment")
doc <- body_add_flextable(doc, two_col_ft(table_results$analysis, "", "Value"), align = "left")
doc <- note(doc, "HR: surgery within 6 months versus no surgery within 6 months. ψ(θ): true CCW log HR (probability limit of the reference analysis) when the Weibull prediction model has on-treatment log HR θ. MCSE: Monte Carlo standard error.")
doc <- body_add_par(doc, "")
doc <- caption(doc, "Table 3: ", "Results (2): performance of the CCW analysis, the simple adjustment and the nested simulation (estimate (MCSE)); MSE = RMSE²")
doc <- body_add_flextable(doc, perf_ft(table_results$performance), align = "left")
doc <- note(doc, paste0(
  "CCW as published: ", n0(B), " datasets per scenario. Simple adjustment: non-parametric correction, CI ",
  "and test learned from these datasets and applied to ", n0(n_val), " independent datasets per scenario. ",
  "Nested simulation: all adjustment steps repeated in each outer dataset (", ns$n_outer, " outer datasets × ",
  ns$n_inner, " inner datasets). MCSE of the bias includes the Monte Carlo error of the true value."
))

doc <- body_add_break(doc)
doc <- caption(doc, "Figure 4: ", "Bias, RMSE, MSE, 95% CI coverage and type I error of the CCW analysis as published, after the simple non-parametric adjustment and in the nested simulation. Points: estimates; bars: 95% Monte Carlo intervals (±1.96 MCSE); dashed lines: nominal values.")
doc <- add_fig(doc, file.path(fig_dir, "fig_performance.png"), 6.5, 6.5 * 4.4 / 15)
doc <- caption(doc, "Figure 5: ", "Histograms of MC estimates: CCW estimates at the calibrated generating value θ* (left) and with no treatment effect (right), with the true values and the non-parametric rejection limits of the 5% test.")
doc <- add_fig(doc, file.path(fig_dir, "fig_distribution.png"), 6.5, 6.5 * 4.6 / 12)
doc <- caption(doc, "Figure 6: ", "Nested simulation: estimated bias correction (left) and bias-adjusted estimate (right) across outer datasets.")
doc <- add_fig(doc, file.path(fig_dir, "fig_nested.png"), 6.5, 6.5 * 4.6 / 12)

# ---- Additional endpoint: survival at 12 months ---------------------------------------
read_ds <- function(nm, ds_name, st = store) {
  tryCatch(targets::tar_read_raw(paste0(nm, "_", ds_name), store = st), error = function(e) NULL)
}
rd_ds <- "maringe_surv12"
table_rd <- read_ds("table_results", rd_ds)
if (!is.null(table_rd)) {
  adj_rd <- read_ds("adjustment", rd_ds)
  pc_rd <- read_ds("perf_compare", rd_ds)
  pp <- function(x, d = 1) formatC(100 * x, format = "f", digits = d)
  doc <- body_add_break(doc)
  doc <- head2(doc, "4.5b Additional endpoint: survival at 12 months")
  doc <- txt(doc, paste0(
    "As in Maringe et al. [35], the effect is also expressed as the difference in survival at 12 months ",
    "between the strategies (surgery minus no surgery), estimated from the same cloned, censored and ",
    "weighted data as the weighted proportion alive at day ", cfg$datasets[[rd_ds]]$pipeline$surv_time,
    " among clones still uncensored, with a patient-clustered SE. The observed difference was ",
    pp(adj_rd$est_obs), " percentage points (95% CI ", pp(adj_rd$ci_lo_obs), " to ", pp(adj_rd$ci_hi_obs), "); ",
    "calibration gave a generating HR of ", f4(exp(read_ds("calibration", rd_ds)$theta_cal)), ". The bias was ",
    pp(adj_rd$bias, 2), " points (MCSE ", pp(adj_rd$bias_mcse, 2), "); the bias-adjusted difference was ",
    pp(adj_rd$est_adj), " points (non-parametric 95% CI ", pp(adj_rd$ci_lo_np), " to ", pp(adj_rd$ci_hi_np), ")."
  ))
  doc <- caption(doc, "Table 3b: ", "Results for the difference in 12-month survival: observed analysis, calibration and QBA adjustment (percentage points)")
  doc <- body_add_flextable(doc, two_col_ft(table_rd$analysis, "", "Value"), align = "left")
  doc <- body_add_par(doc, "")
  doc <- caption(doc, "Table 3c: ", "Performance for the difference in 12-month survival (estimate (MCSE), percentage points)")
  doc <- body_add_flextable(doc, perf_ft(table_rd$performance), align = "left")
}

# ---- 4.6 Additional analyses ------------------------------------------------------
table_rct <- read_ds("table_results_rct", ds)
table_gg <- read_ds("table_results", ds, "_targets_gengamma")
if (!is.null(table_rct) || !is.null(table_gg)) {
  doc <- head2(doc, "4.6 Additional analyses")
}
if (!is.null(table_gg)) {
  adj_gg <- read_ds("adjustment", ds, "_targets_gengamma")
  dgm_gg <- read_ds("dgm", ds, "_targets_gengamma")
  doc <- txt(doc, paste0(
    "As a sensitivity analysis, a more flexible prediction model was applied: the Weibull baseline was ",
    "replaced by a three-parameter generalized gamma baseline (scale σ = ", f2(dgm_gg$sigma), ", shape Q = ",
    f2(dgm_gg$Q), "; Q = 1 is the Weibull; log-likelihood ", f2(dgm_gg$loglik), " versus ",
    f2(dgm_gg$loglik_weibull), "), keeping the proportional-hazards form and the time-varying surgery effect. ",
    "The calibrated generating HR was ", f4(exp(read_ds("calibration", ds, "_targets_gengamma")$theta_cal)),
    " and the bias-adjusted HR ", f4(exp(adj_gg$est_adj)), " (non-parametric 95% CI ", f4(exp(adj_gg$ci_lo_np)),
    "; ", f4(exp(adj_gg$ci_hi_np)), "), against ", f4(exp(adj$est_adj)), " (", f4(exp(adj$ci_lo_np)), "; ",
    f4(exp(adj$ci_hi_np)), ") with the Weibull model. Full results are in Table S3."
  ))
}
if (!is.null(table_rct)) {
  adj_rct <- targets::tar_read_raw(paste0("adjustment_rct_", ds), store = store)
  ns_rct <- targets::tar_read_raw(paste0("nested_rct_", ds), store = store)$summary
  tt_rct <- targets::tar_read_raw(paste0("truth_table_rct_", ds), store = store)
  doc <- txt(doc, paste0(
    "A further analysis does not assume that the analysis pipeline is asymptotically unbiased, which ",
    "requires another way to derive the target estimand: a hypothetical randomised trial in the simulated ",
    "world, analysed without weights. Each patient of the data enters both arms; in the surgery arm surgery ",
    "takes place at the patient's observed date (within the grace period) or, for patients without surgery ",
    "within the grace period, at a date drawn from the observed dates; the arm without surgery is censored ",
    "at day ", cfg$rct$control_censor_day, ". Outcomes are predicted by the prediction model, with a new ",
    "prediction for every copy of a patient in datasets of about ", n_big(cfg$truth$n_target), " patients. ",
    "The trial HR at θ* is ", f4(exp(adj_rct$truth)), " (CCW pipeline limit ", f4(exp(adj$truth)), "), so the ",
    "bias becomes ", f3(adj_rct$bias), " (MCSE ", f3(adj_rct$bias_mcse), "). The adjusted HR is ",
    f4(exp(adj_rct$est_adj)), " (non-parametric 95% CI ", f4(exp(adj_rct$ci_lo_np)), "; ", f4(exp(adj_rct$ci_hi_np)),
    "; nested CI ", f4(exp(ns_rct$ci_lo_np)), "; ", f4(exp(ns_rct$ci_hi_np)), "). These results matter because the ",
    "pipeline remained centred at HR = ", f4(exp(tn$psi_large_sample)), " instead of 1 for large samples under the null (Table X)."
  ))
  doc <- caption(doc, "Table X: ", "Results when not assuming a consistent estimator (target: unweighted hypothetical trial)")
  doc <- body_add_flextable(doc, two_col_ft(table_rct$analysis, "", "Value"), align = "left")
  doc <- body_add_par(doc, "")
  doc <- body_add_flextable(doc, perf_ft(table_rct$performance), align = "left")
}

# ---- Appendix: outer / inner split ---------------------------------------------------
alloc <- read_ds("allocation", ds)
if (!is.null(alloc)) {
  doc <- body_add_break(doc)
  doc <- head2(doc, "Appendix: split of the nested simulation into outer and inner runs (pilot study)")
  doc <- txt(doc, paste0(
    "Pilot: ", alloc$n_outer_pilot, " outer datasets × ", alloc$R_pilot, " inner runs. Spread of the bias-adjusted ",
    "estimate across outer datasets ", f4(alloc$var_obs), " (variance); inner single-run variance σ_w² = ",
    f4(alloc$sigma_w2), "; V̂ = ", f4(alloc$V), "; k = σ_w²/V̂ = ", f2(alloc$k), ". Measured cost: c_o = ",
    f2(alloc$c_o_sec), " s per outer dataset without inner runs, c_i = ", f3(alloc$c_i_sec), " s per inner run ",
    "(including all calibration steps; ", f2(alloc$mean_eval), " steps on average); c_o/c_i = ",
    f1(alloc$cost_ratio), ". The optimum R* = ", f1(alloc$R_star), " inner runs; ", alloc$B_precision,
    " outer datasets give a 2% relative MCSE of the empirical SE, and ", alloc$B_rec, " are used so that the percentile ",
    "limits are stable. Recommended: ", alloc$B_rec, " outer × ", alloc$R_rec, " inner runs (",
    f1(alloc$cost_rec_cpu_h), " CPU hours) instead of ", cfg$nested$outer, " × ", cfg$nested$inner, " (",
    f1(alloc$cost_current_cpu_h), " CPU hours)."
  ))
}

# Supplement
doc <- body_add_break(doc)
doc <- head2(doc, "Online Supplement (draft)")
doc <- caption(doc, "Table S1: ", "Weibull prediction model for death (hazard ratios)")
doc <- body_add_flextable(doc, dgm_ft(table_results$dgm), align = "left")
doc <- note(doc, paste0(dgm$n_events, " deaths. Surgery enters as a time-varying covariate (0 before, 1 from the surgery date)."))
doc <- caption(doc, "Table S2: ", "Coverage of the 95% CI variants")
doc <- body_add_flextable(doc, ci_var_ft(table_results$ci_variants), align = "left")
doc <- note(doc, "Non-parametric: centred MC quantiles added to the bias-adjusted estimate (Section 2.7.2). Reflected: [estimate − q97.5; estimate − q2.5]. Parametric: robust SE multiplied by the CI inflation factor. Below/above: CIs entirely below/above the true value (2.5% each expected).")

if (!is.null(table_gg)) {
  doc <- caption(doc, "Table S3: ", "Sensitivity analysis with a generalized gamma prediction model")
  doc <- body_add_flextable(doc, two_col_ft(table_gg$analysis, "", "Value"), align = "left")
  doc <- body_add_par(doc, "")
  doc <- body_add_flextable(doc, perf_ft(table_gg$performance), align = "left")
}

bw <- balance[weighting == "Weighted"]
bu <- balance[weighting == "Unweighted"]
doc <- caption(doc, "Figure S1: ", paste0(
  "Standardised differences between the surgery and no-surgery strategies among clones at risk at the end ",
  "of the grace period (", bu$n1[1], " and ", bu$n0[1], " clones), before and after censoring weights. ",
  "Dashed lines: ±10%."
))
doc <- add_fig(doc, file.path(fig_dir, "fig_balance.png"), 5.2, 5.2 * 5 / 7)
doc <- note(doc, paste0(
  "Before weighting ", sum(abs(bu$smd) > 10), " of ", nrow(bu), " covariate terms differ by more than 10%; ",
  "after weighting ", sum(abs(bw$smd) > 10), " (largest ", f1(max(abs(bw$smd))), "%). [Maringe et al. report ",
  "all weighted differences within 10% in their registry data; the 200-patient example does not reach this.]"
))

print(doc, target = out_file)
message("Written: ", out_file)
