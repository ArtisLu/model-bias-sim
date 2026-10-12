# ADEMP protocol — modelling-bias QBA, CCW case studies (v0.2)

**Owner:** Artis Luguzis · **Date:** 5 Oct 2026 · **Status:** demo run; to be confirmed with Gerd before the production run.

## A — Aims

1. Apply the draft paper's modelling-bias QBA to the published CCW analysis of two public datasets.
2. Per dataset: does the CCW analysis show bias in the log-HR, under-coverage of its 95% CI, or inflated type I error, when the true effect equals the observed estimate?
3. Compare the two datasets, and show the correction step (paper Section 2.7).

## Base data

| | Maringe et al. (2020, IJE 49:1719) | Webster-Clark et al. (2025, PDS 34:e70223) |
|---|---|---|
| Source | `ccwr::lungcancer` (simulated) | `baseccwcohort.sas7bdat`, CCW_tutorial repository (SynPUF-based) |
| n / events | 200 / 48 deaths within 1 year | 1,590 / 578 readmissions or deaths within 180 days |
| Treatment | Surgery (106 within 6 months; 42 on day 1) | Treatment start after MI discharge (`RX_start`) |
| Covariates | age, sex, performance status, stage, deprivation, Charlson, emergency | age, sex, end-stage renal disease, race (DGM only) |
| Censoring | Administrative at 365.24 days | Truncated at 180 days; early loss to follow-up (31 patients) kept |

Maringe: classic CCW teaching example, already cited in the draft; small sample. SynPUF: the tutorial cohort ships with the repository, so no cohort build from CMS files is needed; larger sample and a different CCW design.

## D — Data-generating mechanism (per dataset)

- **Kept fixed (paper Step 2):** all covariates, each patient's observed treatment date `s` (treatment plan; `s = ∞` if none), and each patient's censoring time.
- **Simulated (Step 3):** event time from a Weibull proportional-hazards model fitted to the original data, with treatment as a **time-varying** exposure: h(t) = λρt^(ρ−1) exp(γ'X + θ·1{t ≥ s}). Treatment is observed only if `s` ≤ event/censoring time.
- Event times by inversion of the cumulative hazard; common random numbers (replicate *b* uses the same uniforms in both scenarios).
- **Scenarios:** θ\* (calibrated so the true CCW effect equals the observed estimate) and θ = 0 (no effect, type I error).

*Why time-varying:* the treatment effect can only start once treatment happens; a DGM with a strategy effect from day 0 would violate "no anticipation".

## E — Estimand and true value

- Estimand: log-HR between the two CCW strategies (per-protocol HR).
  - Maringe: surgery within 6 months vs no surgery within 6 months.
  - SynPUF: treatment start within 90 days vs within 30 days (tutorial's main contrast; the tutorial reports a 180-day risk difference, an HR is used here so both datasets share one estimand type).
- The strategy contrast is not θ (grace-period dilution, non-proportional hazards, non-collapsibility). The true value ψ(θ) is the large-sample limit of the CCW pipeline in the simulated world: the pipeline applied to 24 independent datasets of about 50,000 patients (stacked copies of the base data), averaged.
- **Calibration (paper Step 3):** secant search for θ\* with ψ(θ\*) = observed CCW estimate (12 large datasets per step, tolerance 0.005, |θ| ≤ 3), started at the fitted θ̂.
- **Null scenario (θ = 0):** every patient's outcome is identical under both strategies, so all estimands are exactly 0; the truth is 0. ψ(0) is computed as a diagnostic only: ψ(0) ≠ 0 measures the CCW offset that persists with unlimited data. In the observed-effect scenario the same offset is absorbed into ψ(θ\*), so the bias reported there is the finite-sample part only.

## M — Methods under test (one function `run_pipeline()`, design set per dataset)

**Maringe (`grace_vs_none`).** Clone into both strategies; grace period 182.62 days. Control clones censored at surgery; surgery clones censored at 182.62 days if alive and not operated (deaths before 182.62 without surgery remain in both arms). Control-arm IPCW: Cox model for time to surgery on 7 covariates, weight 1/Ŝc(t|X) updated at each censoring time. Surgery-arm IPCW: logistic model for surgery by 182.62 days among those alive, weight 1/p̂ from 182.62 onwards.

**SynPUF (`window_vs_window`, tutorial code 01–02).** For each window w ∈ {30, 90}: clones followed beyond w without having started are censored at w. Weight 1 in [0, w]. After w: logistic model, among clones censored at w and clones that started in the last 7 days of the window, for being uncensored on age, sex and renal disease; starters in the last 7 days get 1/p̂, earlier starters get 1.

**Both.** Weighted Cox model with strategy only, Efron ties, robust (patient-clustered) SE; Wald 95% CI.

## P — Performance measures

Per dataset and scenario: bias ± MCSE (including the truth's MC error), empirical SE, model SE and relative error, RMSE, 95% CI coverage, rejection rate of H0: HR = 1 (and its direction), failures/warnings. **B = 5,000** replicates per scenario (bias MCSE ≈ 0.004, coverage MCSE ≈ 0.3%). Failures are recorded and excluded; their count is reported. Production target: 10,000.

## Correction (paper Step 7)

- Bias: corrected estimate = observed − bias (the paper's grid search to first order).
- CI: widened by empirical SE / model SE in the observed-effect scenario (when > 1).
- Type I error: p-value = share of null-scenario replicates with |z| ≥ observed |z|; calibrated critical value = 95th percentile of null |z|.

## Change log

| Date | Change | Reason |
|---|---|---|
| 2026-10-02 | v0.1 demo design (Maringe only, θ grid) | First run for colleague demo |
| 2026-10-02 | Control-arm censoring model: weekly pooled logistic → Cox (Maringe) | Weekly discretisation left within-week informative censoring uncorrected |
| 2026-10-02 | Added null-specific measures (coverage of 0, direction of rejections) | Clarify type I error definition |
| 2026-10-05 | v0.2: added SynPUF dataset; θ grid replaced by observed-effect + null scenarios; secant calibration; truth with 24 large datasets; correction = observed − bias | Focus on the observed effect; compare datasets |
