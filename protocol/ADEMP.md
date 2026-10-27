# ADEMP protocol — modelling-bias QBA, CCW case study (v0.5)

**Owner:** Artis Luguzis · **Date:** 27 Oct 2026 · **Status:** paper run (Maringe only, two endpoints), adjustment and additional analyses per manuscript v0.65 (24 Oct 2026). v0.2 (Maringe + SynPUF) is frozen on branch `archive/two-datasets`; SynPUF stays in the code and can be rerun with the `two_datasets` profile.

## A — Aims

1. Apply the draft paper's modelling-bias QBA to the published CCW analysis of the Maringe et al. (2020) lung cancer example.
2. Does the CCW analysis show bias in the log-HR, under-coverage of its 95% CI, or inflated type I error, when the true effect equals the observed estimate?
3. Apply the adjustment step of paper Section 2.7 non-parametrically (inverse calibration and bias correction, CI, test), check it in independent simulated studies, and evaluate the complete adjustment by a nested (double) simulation; present both.
4. (v0.2, `two_datasets` profile) Compare with the SynPUF MI cohort.
5. (v0.5, paper 4.6) Additional analyses: a second endpoint (difference in survival at 12 months); sensitivity of the conclusions to the prediction model (3-parameter generalized gamma instead of Weibull); bias correction towards a hypothetical-trial target for an estimator that is not consistent.
6. (v0.5, paper Appendix) Choose the numbers of outer and inner runs of the nested simulation from a pilot study.

## Base data

| | Maringe et al. (2020, IJE 49:1719) | Webster-Clark et al. (2025, PDS 34:e70223) |
|---|---|---|
| Source | `ccwr::lungcancer` (simulated) | `baseccwcohort.sas7bdat`, CCW_tutorial repository (SynPUF-based) |
| n / events | 200 / 48 deaths within 1 year | 1,590 / 578 readmissions or deaths within 180 days |
| Treatment | Surgery (106 within 6 months; 42 on day 1) | Treatment start after MI discharge (`RX_start`) |
| Covariates | age, sex, performance status, stage, deprivation, Charlson, emergency | age, sex, end-stage renal disease, race (DGM only) |
| Censoring | Administrative at 365.24 days | Administrative at 180 days; loss to follow-up before day 180 in 30 patients |

Maringe: classic CCW teaching example, already cited in the draft; small sample. SynPUF: the tutorial cohort ships with the repository, so no cohort build from CMS files is needed; larger sample and a different CCW design.

## D — Data-generating mechanism (per dataset)

- **Kept fixed (paper Step 2):** all covariates and each patient's observed treatment date `s` (treatment plan; `s = ∞` if none).
- **Simulated (Step 3):** event time from a Weibull proportional-hazards model fitted to the original data, with treatment as a **time-varying** exposure: h(t) = λρt^(ρ−1) exp(γ'X + θ·1{t ≥ s}). Treatment is observed only if `s` ≤ event/censoring time.
- **Censoring (paper Section 2.3; v0.3):** a censoring model is fitted only where there is censoring to model. In Maringe every survivor is followed to 365.24 days (152 of 152), so censoring is administrative only and identical in every replicate (`censoring: none`); a Weibull model with the switched indicator would have no events. For SynPUF (30 losses before day 180) loss to follow-up is simulated from a second Weibull PH model fitted with the event/censoring indicator switched (loss = event; outcome events and day 180 = censored), same covariates, treatment as a time-varying term, so censoring after treatment start has its own hazard. `censoring: fixed` reproduces v0.2 (observed censoring times kept).
- **Sensitivity model (v0.5, paper 4.6):** `dgm: family: gengamma` replaces the Weibull baseline by a generalized gamma with proportional hazards in the covariates (flexsurv parametrisation, μ fixed at 0 with a free intercept, scale σ and shape Q; Weibull at Q = 1), time-varying treatment, fitted by maximum likelihood on time divided by the largest administrative time (conditioning); event times by inversion. Own store (`R_CONFIG_ACTIVE=gengamma`).
- Event and censoring times by inversion of the cumulative hazards; common random numbers (replicate *b* uses the same uniforms in both scenarios; censoring uniforms are drawn after the event uniforms, so event draws are unchanged from v0.2).
- **Scenarios:** θ\* (inverse calibration: the mean CCW estimate over simulated studies of the original size equals the observed estimate) and θ = 0 (no effect, type I error).

*Why time-varying:* the treatment effect can only start once treatment happens; a DGM with a strategy effect from day 0 would violate "no anticipation".

## E — Estimand and true value

- Estimand: log-HR between the two CCW strategies (per-protocol HR). Second endpoint (v0.5, `maringe_surv12`): difference in survival at 12 months between the strategies (arm 1 − arm 0), the weighted proportion alive at day 365 among clones still uncensored (Hájek form), SE from the patient-clustered influence function (weights treated as known); a larger difference favours arm 1, so the "favouring" rates use the opposite tail to the HR.
  - Maringe: surgery within 6 months vs no surgery within 6 months.
  - SynPUF: treatment start within 90 days vs within 30 days (tutorial's main contrast; the tutorial reports a 180-day risk difference, an HR is used here so both datasets share one estimand type).
- The strategy contrast is not θ (grace-period dilution, non-proportional hazards, non-collapsibility). The true value ψ(θ) is the large-sample limit of the CCW pipeline in the simulated world: the pipeline applied to 24 independent datasets of about 50,000 patients (stacked copies of the base data), averaged.
- **Calibration (paper 2.3 and 2.7.1, v0.4):** secant search for θ\* with m(θ\*) = observed CCW estimate, where m(θ) is the mean CCW estimate over 2,000 simulated studies of 200 patients (common random numbers, replicates 1–2,000 of the main run; tolerance 0.002, |θ| ≤ 3), started at θ = −0.7 (v0.5; at the fitted θ̂ before), stopping at once if m(start) is within the tolerance, with a first Newton step using `calib_slope` and secant steps afterwards; tolerances for other endpoints are scaled by `tol_scale`. ψ(θ\*) is then computed as above and the bias b(θ\*) = m(θ\*) − ψ(θ\*) from all B studies at θ\*. (v0.2–v0.3: θ\* with ψ(θ\*) = observed.)
- **Hypothetical-trial target (v0.5, paper 4.6, HR endpoint only):** when the pipeline is not consistent, ψ(θ) is replaced by the log HR of a hypothetical randomised trial in the simulated world, analysed without weights. Implementation assumptions (`rct` in the config): every patient is in both arms; arm 1 gets surgery at the patient's own date if it falls within the grace period, otherwise at a date drawn from the observed in-window dates; arm 0 has no surgery and is censored at day 181 (`control_censor_day`; null = follow-up end); large-sample value from stacked copies as above. Only the truth changes, so the same replicates and inner runs are reused (`R/rct.R`). Check: under the null the target is exactly 0.
- **Null scenario (θ = 0):** every patient's outcome is identical under both strategies, so all estimands are exactly 0; the truth is 0 and no calibration is needed (paper 2.7.3). ψ(0) is computed as a diagnostic only: ψ(0) ≠ 0 measures the CCW offset that persists with unlimited data. At θ\* the same offset is absorbed into ψ(θ\*) (ψ = probability limit of the reference analysis, paper 2.7.1), so the bias reported there is the finite-sample part only.

## M — Methods under test (one function `run_pipeline()`, design set per dataset)

**Maringe (`grace_vs_none`).** Clone into both strategies; grace period 180 days (v0.5; 182.62 before). Control clones censored at surgery; surgery clones censored at day 180 if alive and not operated (deaths before day 180 without surgery remain in both arms). Control-arm IPCW: Cox model for time to surgery on 7 covariates, weight 1/Ŝc(t|X) updated at each censoring time. Surgery-arm IPCW: logistic model for surgery by day 180 among those alive, weight 1/p̂ from day 180 onwards.

**SynPUF (`window_vs_window`, tutorial code 01–02).** For each window w ∈ {30, 90}: clones followed beyond w without having started are censored at w. Weight 1 in [0, w]. After w: logistic model, among clones censored at w and clones that started in the last 7 days of the window, for being uncensored on age, sex and renal disease; starters in the last 7 days get 1/p̂, earlier starters get 1.

**Both (HR endpoint).** Weighted Cox model with strategy only, Efron ties, robust (patient-clustered) SE; Wald 95% CI. **12-month survival (`endpoint: surv12`, Maringe).** Same clones and weights; difference of the weighted survival proportions at day 365, Wald 95% CI.

## P — Performance measures

Per scenario: bias ± MCSE (including the truth's MC error), empirical SE, model SE and relative error, RMSE and MSE = RMSE² (each with MCSE; v0.5), 95% CI coverage and one-sided misses, type I error (rejection rate of H0 at θ = 0, and its direction); power is no longer calculated (v0.5); failures/warnings. **B = 10,000** replicates per scenario (bias MCSE ≈ 0.003, coverage MCSE ≈ 0.2%), plus **10,000 independent replicates** per scenario (separate seed stream) to check the adjusted analysis. Failures are recorded and excluded; their count is reported. The same measures are reported for the simple adjustment (independent replicates) and the nested simulation (outer datasets).

## Adjustment (paper Section 2.7; v0.4)

All adjustments are non-parametric (primary); parametric versions are reported as a secondary check.

- **Bias (2.7.1):** adjusted estimate ψ(θ\*) = ψ̂_obs − b(θ\*), b(θ\*) = m(θ\*) − ψ(θ\*) from the B studies at θ\*; MCSE = √(s²/B + MCSE(ψ)²).
- **CI (2.7.2, non-parametric):** centre the B estimates at θ\* at their mean; CI = adjusted estimate + [q2.5, q97.5] of the centred errors (type-8 quantiles). The reflected form [− q97.5, − q2.5] is reported alongside (exact if the errors are skewed). Parametric: adjusted estimate ± 1.96 · k · SE, k = empirical SE / RMS model SE.
- **Test (2.7.3, non-parametric):** centre the B no-effect estimates at their mean m0; reject H0 if ψ̂_obs − m0 is below q2.5 or above q97.5 of the centred no-effect estimates (separate lower and upper limits); simulation-calibrated p = min{1, 2 · min(P(e0 ≤ d), P(e0 ≥ d))}. Parametric: centred |z| with critical value = 95th percentile under θ = 0.
- **Check (2.7.2 step 5c, 2.7.3 step 9):** the adjustment learned from the main replicates is applied unchanged to 10,000 independent replicates per scenario; bias, RMSE, coverage (and misses by side), type I error and its direction are reported with MCSE. Replaces the two-fold cross-fitting of v0.3.
- **Nested simulation (2.7.2):** outer datasets = the first 1,250 independent replicates (v0.5; 500 before) (at θ\* and at θ = 0). Per outer dataset at θ\*: CCW estimate; refit of the Weibull prediction model; inverse calibration of θ_r with 50 inner studies per step (v0.5; 250 before) (Newton step from θ\* with the slope of the main calibration, then secant; common random numbers within the outer dataset); ψ_r(θ_r) under the refitted model (2 large datasets of 25,000); b_r = m_r − ψ_r; adjusted estimate ψ̂_r − b_r; non-parametric CI from the inner studies at θ_r. Coverage = share of outer CIs containing ψ(θ\*). Per outer dataset at θ = 0: refit, 50 inner no-effect studies, own rejection limits; type I error = share of rejections. Inner seeds are derived from the outer replicate's auxiliary seed. Reported: coverage, type I error, bias/empirical SE/RMSE of the outer adjusted estimates, SD of b_r and its within-dataset MC error, failures, calibrations at the bound, computing time; paired comparison with the simple adjustment on the same outer datasets. **Nested CI of the observed analysis:** adjusted estimate + centred 2.5/97.5% quantiles of the outer adjusted estimates (includes the variation of the estimated bias correction).
- **Outer / inner split (paper Appendix, v0.5; `R/allocation.R`):** a pilot (`pilot:` in the config, 50 outer × 50 inner by default) gives the ANOVA decomposition Var_obs = V + σ_w²/R of the bias-adjusted estimates (σ_w² within outer datasets, V without inner noise), k = σ_w²/V and the costs c_o (outer replicate without inner runs, incl. the large-sample truth) and c_i (one inner run, all calibration steps). R\* = (k/2)(1 + √(1 + 8c_o/(k c_i))) ≈ √(2k c_o/c_i); B for a relative MCSE of the empirical SE of 2% is 1 + 1/(2·0.02²) ≈ 1,250; for percentile limits at least 1,000 outer runs and R, B scaled with the budget as C^(1/3) and C^(2/3). Check that V = Var_obs(r) − σ_w²/r does not depend on r (first r inner runs). `nested: shrink: true` shrinks the outer deviations by √(V/Var_obs). **Pilot of 27 Oct 2026 (Maringe, log HR; 50 × 50, seed_pilot 20261027):** Var_obs = 0.0965, σ_w² = 0.0771, V = 0.0950, k = 0.81; c_o = 8.8 s (refit, truth of 2 × 25,000, outer calibration overhead), c_i = 0.126 s (one inner study including the 2.5 calibration steps on average), c_o / c_i = 70; R\* = 11.1; B = 1,251 for a 2% relative MCSE of the empirical SE; shrink factor 0.992. V is stable across r = 5…50 (0.090–0.095), as the decomposition requires. **Production: 1,250 outer × 50 inner**, not 1,250 × 11: R\* minimises the variance of V̂, but the limits of each outer CI and test are quantiles of its inner studies, and 11 of them are 25% too narrow (mean width 0.85 against 1.13 log HR; 25 inner: 1.11), so coverage and type I error of the outer CIs would be biased. 1,250 × 50 costs 5.3 CPU h (500 × 250: 5.6 CPU h) and has a 2.6 times lower variance of V̂ at equal cost. Profile `nested_rstar` runs 1,250 × 11 (bias correction, SE and RMSE of the outer estimates only).
- **Optional (config/belt.yml, off by default since v0.4):** simulation bootstrap-t and confidence belt.

## Additional analyses (v0.5, paper 4.6)

| Analysis | Config | Output |
|---|---|---|
| Survival difference at 12 months | `maringe_surv12` in `active_datasets` | Same tables and figures as the HR endpoint, in percentage points (MSE in squared points) |
| Generalized-gamma prediction model | profile `gengamma` | Own store; HR and survival difference, compared with the Weibull run |
| Hypothetical-trial target | `rct:` | `truth_rct`, `adjustment_rct`, `nested_rct`, `perf_compare_rct` (HR only) |
| Outer / inner split | `pilot:` | `allocation`, `pilot_anova` |

## Change log

| Date | Change | Reason |
|---|---|---|
| 2026-10-02 | v0.1 demo design (Maringe only, θ grid) | First run for colleague demo |
| 2026-10-02 | Control-arm censoring model: weekly pooled logistic → Cox (Maringe) | Weekly discretisation left within-week informative censoring uncorrected |
| 2026-10-02 | Added null-specific measures (coverage of 0, direction of rejections) | Clarify type I error definition |
| 2026-10-05 | v0.2: added SynPUF dataset; θ grid replaced by observed-effect + null scenarios; secant calibration; truth with 24 large datasets; correction = observed − bias | Focus on the observed effect; compare datasets |
| 2026-10-12 | v0.3: Maringe only for the paper (SynPUF kept in code, `two_datasets` profile); Weibull censoring model with switched indicator and time-varying treatment where loss to follow-up exists (SynPUF; Maringe has administrative censoring only); B = 10,000; cross-fitted performance of the adjusted analysis; paper tables/figures | Paper uses one dataset; paper Section 2.3 (censoring, including after treatment start, is predicted) |
| 2026-10-13 | Added simulation bootstrap-t CI as an alternative; CI comparison table and figure | Wald misses are lopsided (est–SE correlation), which symmetric widening cannot fix |
| 2026-10-13 | Added confidence belt (test inversion over a θ grid) as the nominal-coverage CI; settings in config/belt.yml | Bootstrap-t assumes T is pivotal; its quantiles drift with the true effect |
| 2026-10-22 | v0.4: inverse calibration m(θ\*) = observed; adjusted estimate ψ(θ\*) = observed − b(θ\*); non-parametric CI and test (centred MC quantiles), parametric as secondary; independent validation replicates (seed_validation) instead of cross-fitting; nested (double) simulation for CI and type I error; bootstrap-t and belt made optional | Manuscript v0.64 (17 Oct 2026) Section 2.7 and the instructions of 22 Oct 2026 |
| 2026-10-22 | Decisions confirmed: (1) truth ψ(θ) = probability limit of the CCW analysis (paper 2.7.1 formula), not a counterfactual strategy-population contrast, so the large-sample CCW offset ψ(0) ≠ 0 is not counted as bias; (2) the non-parametric CI adds the centred quantiles (paper form) and is the primary CI, with the reflected form reported alongside as a check | Kept as implemented in v0.4 |
| 2026-10-27 | v0.5: calibration starts at θ = −0.7; grace period exactly 180 days; RMSE and MSE in all outputs, power removed; second endpoint (12-month survival difference); generalized-gamma sensitivity profile; hypothetical-trial target for a non-consistent estimator (assumptions listed in E); pilot study and allocation of outer / inner runs (paper Appendix), optional shrinkage of outer deviations | Manuscript v0.65 (24 Oct 2026) Sections 4.6 and Appendix and the instructions of 27 Oct 2026 |
