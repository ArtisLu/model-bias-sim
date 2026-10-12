# Word draft of paper Section 4 (Methods text, Results text, tables, figures) from the
# targets store. Run from the project root after targets::tar_make():
#   source("reports/paper/build_section4.R")
# Kept outside the pipeline while the wording is being iterated; reads run times from tar_meta().

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
for (nm in c("table_baseline", "table_results", "perf_adjusted", "performance", "correction",
             "calibration", "dgm", "obs_result", "base_data", "truth_table", "sim_results",
             "balance")) {
  assign(nm, targets::tar_read_raw(paste0(nm, "_", ds), store = store))
}
spec <- cfg$datasets[[ds]]
fig_dir <- file.path(out_root, "paper", ds)
out_file <- file.path(out_root, "paper", "section4_draft.docx")

# ---- Numbers used in the text --------------------------------------------------
f1 <- function(x) formatC(x, format = "f", digits = 1)
f2 <- function(x) formatC(x, format = "f", digits = 2)
f3 <- function(x) formatC(x, format = "f", digits = 3)
pc <- function(x, d = 1) paste0(formatC(100 * x, format = "f", digits = d), "%")
n0 <- function(x) format(x, big.mark = ",")
hr_ci <- function(e, lo, hi) paste0(f2(exp(e)), " (95% CI ", f2(exp(lo)), "–", f2(exp(hi)), ")")
po <- performance[scenario == "observed"]
pn <- performance[scenario == "null"]
cr <- correction
cal <- calibration
bd <- base_data
tn <- truth_table[scenario == "null"]
adj <- function(m) perf_adjusted[analysis == "CCW + QBA adjustment" & measure == m]
mc_ci <- function(v, mc) paste0(f3(v - 1.96 * mc), " to ", f3(v + 1.96 * mc))
B <- po$B
n_fail <- po$n_fail + pn$n_fail
n_trt <- sum(is.finite(bd$trt_time))
# Patterns repeat the time of their branches, so only stems and branches are summed
meta <- targets::tar_meta(fields = c(name, type, time, seconds), targets_only = TRUE, store = store)
meta <- meta[meta$type %in% c("stem", "branch"), ]
secs_sim <- sum(meta[grepl(paste0("^sim_results_", ds, "_"), meta$name), ]$seconds, na.rm = TRUE)
secs_all <- sum(meta$seconds, na.rm = TRUE)
heavy <- meta[grepl(paste0("^(calibration|truth|sim_results)_", ds), meta$name), ]
mins_wall <- as.numeric(difftime(max(heavy$time), min(heavy$time - heavy$seconds), units = "mins"))
rs <- targets::tar_read(run_settings, store = store)
r_ver <- paste(R.version$major, R.version$minor, sep = ".")
n_surg_day1 <- sum(bd$trt_time <= 1)

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
                      adjusted = "CCW + QBA adjustment") |>
    align(j = 2:3, align = "center", part = "all") |>
    style_table() |>
    width(j = 1:3, width = c(3.5, 1.6, 1.6))
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
doc <- head1(doc, "Section 4 – Simulation example (draft for insertion)")
doc <- note(doc, paste0(
  "Generated by reports/paper/build_section4.R from the results store ", store, " (run of ",
  format(Sys.Date(), "%d %b %Y"), ", B = ", n0(B), " per scenario, protocol v0.3). ",
  "Text in square brackets needs a decision or a check by the authors. Numbering of tables ",
  "and figures follows the current draft (Table 1–2 and Figure 1 are already used)."
))

# 4.1 Used dataset
doc <- head2(doc, "4.1 Used dataset")
doc <- txt(doc, paste0(
  "We reanalyse the example dataset of Maringe et al. [36] to check whether the applied clone-censor-",
  "weighting (CCW) approach is unbiased, has nominal 95% CI coverage and no type I error inflation. ",
  "The dataset comprises ", n0(nrow(bd)), " simulated elderly patients with stage I or II lung cancer, ",
  "generated by the authors to mimic their cancer registry study and distributed with the tutorial ",
  "[it is the simulated supplementary dataset, not the registry data of their Table S2; data from the ccwr ",
  "R package]. Time zero is the date of diagnosis and patients are followed for one year. The ",
  "treatment of interest is major surgery, received by ", n0(n_trt), " patients (", pc(n_trt / nrow(bd)),
  "), ", n_surg_day1, " of them on day 1; the median time to surgery among operated patients was ",
  stats::median(bd$trt_time[is.finite(bd$trt_time)]), " days. The endpoint is death within one ",
  "year, observed in ", n0(sum(bd$status)), " patients (", pc(mean(bd$status)), "). All patients ",
  "alive at one year were followed until day 365, so the only censoring is administrative. Baseline ",
  "covariates are age, sex, performance status, TNM stage, deprivation score, Charlson comorbidity ",
  "index and emergency presentation (Table 3)."
))
doc <- caption(doc, "Table 3: ", "Patient characteristics of the Maringe example dataset by surgery during follow-up")
doc <- body_add_flextable(doc, baseline_ft(table_baseline), align = "left")
doc <- note(doc, "Standardised mean differences: one value for binary characteristics, one per category for performance status and Charlson index, continuous scale for age and deprivation score. Deprivation score as coded in the simulated dataset (0–1 scale).")

# 4.2 Model
doc <- head2(doc, "4.2 Model")
doc <- txt(doc, paste0(
  "Reference model. The reference model is the CCW analysis of Maringe et al. [36], comparing the ",
  "strategies 'surgery within 6 months of diagnosis' and 'no surgery within 6 months' (grace period ",
  spec$pipeline$grace, " days). Each patient is cloned into both strategies. Clones in the no-surgery ",
  "strategy are censored on the day of surgery if surgery happens within the grace period. Clones in ",
  "the surgery strategy are censored at the end of the grace period if they are alive and have not been ",
  "operated; patients who die within the grace period without surgery remain in both strategies. First, ",
  "the weights are calculated. In the no-surgery strategy, the probability of remaining uncensored is ",
  "estimated from a Cox model for time to surgery and the weights are updated at every surgery date; in ",
  "the surgery strategy, a logistic model for surgery within the grace period among patients alive at its ",
  "end gives the weights from that time onwards. Both models include age, sex, performance status, stage, ",
  "deprivation score, Charlson index and emergency presentation. The strategies are compared by a ",
  "weighted Cox model with the strategy as the only covariate, Efron ties and a robust, patient-clustered ",
  "variance. Figure 4 shows how example patients are cloned, censored and weighted, and Figure S2 the ",
  "covariate balance between the strategies at 6 months before and after weighting."
))
doc <- txt(doc, paste0(
  "Preserved data. Not only the baseline data are maintained but also each patient's observed date of ",
  "surgery; patients who were not operated keep no surgery date. Surgery is counted as received only if ",
  "its date precedes the simulated end of follow-up."
))
doc <- txt(doc, paste0(
  "Weibull prediction model. Times to death are predicted from a Weibull proportional hazards model ",
  "fitted to the original data by maximum likelihood, with the same seven covariates and surgery as a ",
  "time-varying covariate, h(t | X) = λρt^(ρ−1) exp(γ′X + θ·1{t ≥ s}), where s is the patient's surgery ",
  "date. The effect of surgery can therefore act only after surgery has taken place. Death times are ",
  "drawn by inverting each patient's cumulative hazard; patients still alive at day 365 are censored ",
  "there. Because all survivors in the original data were followed for the full year, there is no loss ",
  "to follow-up to model in this example: the censoring model of Section 2.3 reduces to administrative ",
  "censoring at one year. [For datasets with loss to follow-up, the simulation code fits a second Weibull ",
  "model with the event/censoring indicator switched and treatment as a time-varying covariate.] ",
  "Parameter estimates are given in Table S1 and the fit to the observed data in Figure S1."
))

# 4.3 Calibration
doc <- head2(doc, "4.3 Calibration")
doc <- txt(doc, paste0(
  "The hazard ratio parameter θ in the Weibull prediction model describes the treatment effect of ",
  "surgery once received, at fixed covariates. The CCW hazard ratio, on the other hand, compares ",
  "treatment strategies over the follow-up period and differs from θ because of grace-period dilution ",
  "(in the surgery strategy, patients are untreated until their surgery date, so the strategy hazard ",
  "ratio starts at 1 and moves towards θ over time) and because it is a marginal rather than a ",
  "conditional contrast."
))
doc <- txt(doc, paste0(
  "The true CCW hazard ratio in the simulated world, ψ(θ), was defined as the large-sample limit of the ",
  "CCW analysis: the reference model was applied to ", cfg$truth$n_rep, " independent simulated datasets ",
  "of ", n0(ceiling(cfg$truth$n_target / nrow(bd)) * nrow(bd)), " patients each (stacked copies of the ",
  "baseline data) and the estimates were averaged. The calibration was performed by a secant search ",
  "starting from the fitted coefficient θ̂ = ", f2(dgm$theta_hat), " (HR ", f2(exp(dgm$theta_hat)),
  "): θ was tuned until ψ(θ) equalled the observed CCW hazard ratio of ", f2(exp(obs_result$estimate)),
  ". This needed θ* = ", f2(cal$theta_star), " (Weibull HR ", f2(exp(cal$theta_star)), "), which results ",
  "in HR_CCW = ", f2(exp(po$truth)), " (", nrow(cal$path) - 2, " search step",
  if (nrow(cal$path) - 2 == 1) "" else "s", " after the two starting values; MCSE of the true log HR ",
  f3(po$truth_mcse), "). [Current draft: 'θ = −0.70, which results ",
  "in HR_CCW = 0.48'; HR_CCW equals the observed HR by construction, and the Weibull HR is exp(θ*) = ",
  f2(exp(cal$theta_star)), ".]"
))
doc <- txt(doc, paste0(
  "A calibration step is not needed for the type I error rate, as conditional and marginal treatment ",
  "effects are both HR = 1 in this case. With unlimited data, the CCW analysis would nevertheless ",
  "estimate HR ", f2(exp(tn$psi_large_sample)), " when surgery has no effect, a small offset that the ",
  "finite-sample results below include."
))

# 4.4 Simulation and performance measures
doc <- head2(doc, "4.4 Simulation settings and performance measures")
doc <- txt(doc, paste0(
  "For each of the two scenarios (θ = θ* and θ = 0), ", n0(B), " datasets were simulated and reanalysed ",
  "with the reference model. Replicate b used the same random numbers in both scenarios (common random ",
  "numbers). Performance measures were the bias of the log HR, the empirical and the model-based ",
  "(robust) standard error, the coverage of the nominal 95% CI and the rejection rate of H0: HR = 1 ",
  "(power under θ*, type I error under θ = 0), each with its Monte Carlo standard error (MCSE); the ",
  "MCSE of the bias includes the Monte Carlo error of the true value. Failed model fits were recorded ",
  "and excluded (n = ", n_fail, ")."
))
doc <- txt(doc, paste0(
  "Adjustment (Section 2.7). The bias-corrected estimate is the observed estimate minus the estimated ",
  "bias, the confidence interval is widened by the ratio of empirical to model-based SE when this ratio ",
  "exceeds 1, and the critical value of the test is replaced by the 95th percentile of |z| in the null ",
  "scenario. To check that the adjusted analysis performs as intended without judging it on the ",
  "replicates that defined it, the adjustments were estimated on one half of the replicates and ",
  "applied to the other half (two-fold cross-fitting), in both scenarios."
))
doc <- txt(doc, paste0(
  "Software and run time. The simulation was programmed in R ", r_ver, " (packages survival, data.table ",
  "and targets; parallel workers via crew) and is available in the Online Supplement. The simulation ",
  "took ", f1(mins_wall), " minutes with ", rs$workers, " parallel worker", if (rs$workers > 1) "s" else "",
  " on a cloud virtual machine (", f1(secs_all / 60), " minutes of computing time, of which ",
  f1(secs_sim / 60), " minutes for the 2 × ", n0(B), " replicates, ", f2(secs_sim / (2 * B)),
  " seconds per replicate)."
))

# 4.5 Results
doc <- head2(doc, "4.5 Results")
bias_sig <- abs(po$bias) > 1.96 * po$bias_mcse_total
cov_sig <- abs(po$coverage - 0.95) > 1.96 * po$coverage_mcse
t1_sig <- abs(pn$reject_0 - 0.05) > 1.96 * pn$reject_0_mcse
doc <- txt(doc, paste0(
  "The observed CCW analysis estimated HR = ", hr_ci(obs_result$estimate, obs_result$ci_lo,
  obs_result$ci_hi), " for surgery within 6 months versus no surgery within 6 months (Table 4). With the ",
  "true CCW HR set to this value, the average estimate across ", n0(B), " simulated studies was HR ",
  f2(exp(po$mean_est)), ": a bias of ", f3(po$bias), " on the log scale (95% MC interval ",
  mc_ci(po$bias, po$bias_mcse_total), "; ", f1(abs(100 * po$std_bias)), "% of the empirical SE), ",
  if (bias_sig) "which is statistically detectable but small. " else "which is not statistically different from zero. ",
  "The robust SE was on average ", f1(abs(100 * po$rel_err_modse)), "% ",
  if (po$rel_err_modse < 0) "smaller" else "larger", " than the empirical SE, and the coverage of the ",
  "nominal 95% CI was ", pc(po$coverage), " (MCSE ", pc(po$coverage_mcse), ")",
  if (cov_sig) ", a significant departure from 95%. " else ", consistent with 95%. ",
  "With no effect of surgery, ", pc(pn$reject_0), " (MCSE ", pc(pn$reject_0_mcse), ") of the simulated ",
  "studies rejected H0 at the 5% level", if (t1_sig) ", a significant inflation of the type I error" else "",
  "; ", pc(pn$reject_0_lower), " were rejections in favour of surgery. Power to detect the observed ",
  "effect was ", pc(po$reject_0), " (Table 5, Figures 2 and 3)."
))
doc <- txt(doc, paste0(
  "The QBA-adjusted analysis gave HR = ", hr_ci(cr$est_corr, cr$ci_lo_corr, cr$ci_hi_corr), ", with ",
  "a CI inflation factor of ", f2(cr$se_inflation), " and a calibrated critical value of |z| = ",
  f2(cr$z_crit_calibrated), " instead of 1.96; the simulation-calibrated p-value was ",
  formatC(cr$p_calibrated, format = "f", digits = 3), " (observed p = ",
  formatC(obs_result$p_value, format = "f", digits = 3), "). In the cross-fitted evaluation the ",
  "adjusted analysis had 95% CI coverage of ",
  pc(adj("95% CI coverage")$value), " and a type I error of ", pc(adj("Type I error")$value),
  "; its bias is zero by construction (Table 5, Figure 2)."
))
doc <- note(doc, paste0(
  "[Discussion, 2nd paragraph: the current draft quotes type I error 7.1% and coverage 93.7%. This run ",
  "(B = ", n0(B), "): bias ", f3(po$bias), ", coverage ", pc(po$coverage), ", type I error ",
  pc(pn$reject_0), ".]"
))

doc <- caption(doc, "Table 4: ", "Observed CCW analysis, prediction model and calibration")
doc <- body_add_flextable(doc, two_col_ft(table_results$analysis, "", "Value"), align = "left")
doc <- note(doc, "HR: surgery within 6 months versus no surgery within 6 months. ψ(θ): large-sample CCW log HR when the Weibull prediction model has on-treatment log HR θ. MCSE: Monte Carlo standard error.")
doc <- body_add_par(doc, "")
doc <- caption(doc, "Table 5: ", paste0("Performance of the CCW analysis before and after QBA adjustment (", n0(B), " simulated studies per scenario; estimate (MCSE))"))
doc <- body_add_flextable(doc, perf_ft(table_results$performance), align = "left")
doc <- note(doc, "Adjusted analysis: bias correction, CI widening and calibrated critical value, each estimated on the other half of the replicates (two-fold cross-fitting). MCSE of the bias includes the Monte Carlo error of the true value.")

doc <- body_add_break(doc)
doc <- caption(doc, "Figure 2: ", "Bias, 95% CI coverage, type I error and power of the CCW analysis before and after QBA adjustment. Points: estimates; bars: 95% Monte Carlo intervals (±1.96 MCSE); dashed lines: nominal values. The adjusted bias is zero by construction.")
doc <- add_fig(doc, file.path(fig_dir, "fig_performance.png"), 6.5, 6.5 * 4.2 / 13)
doc <- caption(doc, "Figure 3: ", "Distribution of the simulated CCW estimates when the true CCW HR equals the observed estimate (left), and of the z statistics when there is no treatment effect (right), with nominal and calibrated critical values.")
doc <- add_fig(doc, file.path(fig_dir, "fig_distribution.png"), 6.5, 6.5 * 4.6 / 12)
doc <- caption(doc, "Figure 4: ", "Cloning, artificial censoring and weighting in the CCW analysis, for five example patients. Shaded: grace period.")
doc <- add_fig(doc, file.path(out_root, paste0("fig_protocol_", ds, ".png")), 6.3, 6.3 * 4.6 / 8.5)

# Supplement
doc <- body_add_break(doc)
doc <- head2(doc, "Online Supplement (draft)")
doc <- caption(doc, "Table S1: ", "Weibull prediction model for death (hazard ratios)")
doc <- body_add_flextable(doc, dgm_ft(table_results$dgm), align = "left")
doc <- note(doc, paste0(dgm$n_events, " deaths. Surgery enters as a time-varying covariate (0 before, 1 from the surgery date)."))
doc <- caption(doc, "Figure S1: ", "Fit of the Weibull prediction model: Kaplan–Meier estimate in the original data versus the survival predicted by the Weibull model (averaged over patients, with each patient's observed surgery date).")
doc <- add_fig(doc, file.path(fig_dir, "fig_model_fit.png"), 4.6, 4.6 * 4.6 / 6)

bw <- balance[weighting == "Weighted"]
bu <- balance[weighting == "Unweighted"]
doc <- caption(doc, "Figure S2: ", paste0(
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
