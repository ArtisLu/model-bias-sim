# ADEMP protocol — modelling-bias QBA, CCW case study (v0.4)

**Owner:** Artis Luguzis · **Date:** 22 Oct 2026 · **Status:** paper run (Maringe only), adjustment per manuscript v0.64 (17 Oct 2026). v0.2 (Maringe + SynPUF) is frozen on branch `archive/two-datasets`; SynPUF stays in the code and can be rerun with the `two_datasets` profile.

## A — Aims

1. Apply the draft paper's modelling-bias QBA to the published CCW analysis of the Maringe et al. (2020) lung cancer example.
2. Does the CCW analysis show bias in the log-HR, under-coverage of its 95% CI, or inflated type I error, when the true effect equals the observed estimate?
3. Apply the adjustment step of paper Section 2.7 non-parametrically (inverse calibration and bias correction, CI, test), check it in independent simulated studies, and evaluate the complete adjustment by a nested (double) simulation; present both.
4. (v0.2, `two_datasets` profile) Compare with the SynPUF MI cohort.

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
- Event and censoring times by inversion of the cumulative hazards; common random numbers (replicate *b* uses the same uniforms in both scenarios; censoring uniforms are drawn after the event uniforms, so event draws are unchanged from v0.2).
- **Scenarios:** θ\* (inverse calibration: the mean CCW estimate over simulated studies of the original size equals the observed estimate) and θ = 0 (no effect, type I error).

*Why time-varying:* the treatment effect can only start once treatment happens; a DGM with a strategy effect from day 0 would violate "no anticipation".

## E — Estimand and true value

- Estimand: log-HR between the two CCW strategies (per-protocol HR).
  - Maringe: surgery within 6 months vs no surgery within 6 months.
  - SynPUF: treatment start within 90 days vs within 30 days (tutorial's main contrast; the tutorial reports a 180-day risk difference, an HR is used here so both datasets share one estimand type).
- The strategy contrast is not θ (grace-period dilution, non-proportional hazards, non-collapsibility). The true value ψ(θ) is the large-sample limit of the CCW pipeline in the simulated world: the pipeline applied to 24 independent datasets of about 50,000 patients (stacked copies of the base data), averaged.
- **Calibration (paper 2.3 and 2.7.1, v0.4):** secant search for θ\* with m(θ\*) = observed CCW estimate, where m(θ) is the mean CCW estimate over 2,000 simulated studies of 200 patients (common random numbers, replicates 1–2,000 of the main run; tolerance 0.002, |θ| ≤ 3), started at the fitted θ̂. ψ(θ\*) is then computed as above and the bias b(θ\*) = m(θ\*) − ψ(θ\*) from all B studies at θ\*. (v0.2–v0.3: θ\* with ψ(θ\*) = observed.)
- **Null scenario (θ = 0):** every patient's outcome is identical under both strategies, so all estimands are exactly 0; the truth is 0 and no calibration is needed (paper 2.7.3). ψ(0) is computed as a diagnostic only: ψ(0) ≠ 0 measures the CCW offset that persists with unlimited data. At θ\* the same offset is absorbed into ψ(θ\*) (ψ = probability limit of the reference analysis, paper 2.7.1), so the bias reported there is the finite-sample part only.

## M — Methods under test (one function `run_pipeline()`, design set per dataset)

**Maringe (`grace_vs_none`).** Clone into both strategies; grace period 182.62 days. Control clones censored at surgery; surgery clones censored at 182.62 days if alive and not operated (deaths before 182.62 without surgery remain in both arms). Control-arm IPCW: Cox model for time to surgery on 7 covariates, weight 1/Ŝc(t|X) updated at each censoring time. Surgery-arm IPCW: logistic model for surgery by 182.62 days among those alive, weight 1/p̂ from 182.62 onwards.

**SynPUF (`window_vs_window`, tutorial code 01–02).** For each window w ∈ {30, 90}: clones followed beyond w without having started are censored at w. Weight 1 in [0, w]. After w: logistic model, among clones censored at w and clones that started in the last 7 days of the window, for being uncensored on age, sex and renal disease; starters in the last 7 days get 1/p̂, earlier starters get 1.

**Both.** Weighted Cox model with strategy only, Efron ties, robust (patient-clustered) SE; Wald 95% CI.

## P — Performance measures

Per scenario: bias ± MCSE (including the truth's MC error), empirical SE, model SE and relative error, RMSE, 95% CI coverage and one-sided misses, rejection rate of H0: HR = 1 (power at θ\*, type I error at θ = 0, and its direction), failures/warnings. **B = 10,000** replicates per scenario (bias MCSE ≈ 0.003, coverage MCSE ≈ 0.2%), plus **10,000 independent replicates** per scenario (separate seed stream) to check the adjusted analysis. Failures are recorded and excluded; their count is reported. The same measures are reported for the simple adjustment (independent replicates) and the nested simulation (outer datasets).

## Adjustment (paper Section 2.7; v0.4)

All adjustments are non-parametric (primary); parametric versions are reported as a secondary check.

- **Bias (2.7.1):** adjusted estimate ψ(θ\*) = ψ̂_obs − b(θ\*), b(θ\*) = m(θ\*) − ψ(θ\*) from the B studies at θ\*; MCSE = √(s²/B + MCSE(ψ)²).
- **CI (2.7.2, non-parametric):** centre the B estimates at θ\* at their mean; CI = adjusted estimate + [q2.5, q97.5] of the centred errors (type-8 quantiles). The reflected form [− q97.5, − q2.5] is reported alongside (exact if the errors are skewed). Parametric: adjusted estimate ± 1.96 · k · SE, k = empirical SE / RMS model SE.
- **Test (2.7.3, non-parametric):** centre the B no-effect estimates at their mean m0; reject H0 if ψ̂_obs − m0 is below q2.5 or above q97.5 of the centred no-effect estimates (separate lower and upper limits); simulation-calibrated p = min{1, 2 · min(P(e0 ≤ d), P(e0 ≥ d))}. Parametric: centred |z| with critical value = 95th percentile under θ = 0.
- **Check (2.7.2 step 5c, 2.7.3 step 9):** the adjustment learned from the main replicates is applied unchanged to 10,000 independent replicates per scenario; bias, coverage (and misses by side), power, type I error and its direction are reported with MCSE. Replaces the two-fold cross-fitting of v0.3.
- **Nested simulation (2.7.2):** outer datasets = the first 500 independent replicates (at θ\* and at θ = 0). Per outer dataset at θ\*: CCW estimate; refit of the Weibull prediction model; inverse calibration of θ_r with 250 inner studies per step (Newton step from θ\* with the slope of the main calibration, then secant; common random numbers within the outer dataset); ψ_r(θ_r) under the refitted model (2 large datasets of 25,000); b_r = m_r − ψ_r; adjusted estimate ψ̂_r − b_r; non-parametric CI from the inner studies at θ_r. Coverage = share of outer CIs containing ψ(θ\*). Per outer dataset at θ = 0: refit, 250 inner no-effect studies, own rejection limits; type I error = share of rejections. Inner seeds are derived from the outer replicate's auxiliary seed. Reported: coverage, type I error, bias/empirical SE/RMSE of the outer adjusted estimates, SD of b_r and its within-dataset MC error, failures, calibrations at the bound, computing time; paired comparison with the simple adjustment on the same outer datasets. **Nested CI of the observed analysis:** adjusted estimate + centred 2.5/97.5% quantiles of the outer adjusted estimates (includes the variation of the estimated bias correction).
- **Optional (config/belt.yml, off by default since v0.4):** simulation bootstrap-t and confidence belt.

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
