# ADEMP protocol — modelling-bias QBA, CCW case study (v0.3)

**Owner:** Artis Luguzis · **Date:** 12 Oct 2026 · **Status:** paper run (Maringe only). v0.2 (Maringe + SynPUF) is frozen on branch `archive/two-datasets`; SynPUF stays in the code and can be rerun with the `two_datasets` profile.

## A — Aims

1. Apply the draft paper's modelling-bias QBA to the published CCW analysis of the Maringe et al. (2020) lung cancer example.
2. Does the CCW analysis show bias in the log-HR, under-coverage of its 95% CI, or inflated type I error, when the true effect equals the observed estimate?
3. Show the adjustment step (paper Section 2.7) and check that the adjusted analysis restores nominal performance.
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
- **Scenarios:** θ\* (calibrated so the true CCW effect equals the observed estimate) and θ = 0 (no effect, type I error).

*Why time-varying:* the treatment effect can only start once treatment happens; a DGM with a strategy effect from day 0 would violate "no anticipation".

## E — Estimand and true value

- Estimand: log-HR between the two CCW strategies (per-protocol HR).
  - Maringe: surgery within 6 months vs no surgery within 6 months.
  - SynPUF: treatment start within 90 days vs within 30 days (tutorial's main contrast; the tutorial reports a 180-day risk difference, an HR is used here so both datasets share one estimand type).
- The strategy contrast is not θ (grace-period dilution, non-proportional hazards, non-collapsibility). The true value ψ(θ) is the large-sample limit of the CCW pipeline in the simulated world: the pipeline applied to 24 independent datasets of about 50,000 patients (stacked copies of the base data), averaged.
- **Calibration (paper Step 3):** secant search for θ\* with ψ(θ\*) = observed CCW estimate (12 large datasets per step, tolerance 0.005, |θ| ≤ 3), started at the fitted θ̂ and at 0 (needed anyway for the null scenario); the large datasets of a step run in parallel and are reused for the true values.
- **Null scenario (θ = 0):** every patient's outcome is identical under both strategies, so all estimands are exactly 0; the truth is 0. ψ(0) is computed as a diagnostic only: ψ(0) ≠ 0 measures the CCW offset that persists with unlimited data. In the observed-effect scenario the same offset is absorbed into ψ(θ\*), so the bias reported there is the finite-sample part only.

## M — Methods under test (one function `run_pipeline()`, design set per dataset)

**Maringe (`grace_vs_none`).** Clone into both strategies; grace period 182.62 days. Control clones censored at surgery; surgery clones censored at 182.62 days if alive and not operated (deaths before 182.62 without surgery remain in both arms). Control-arm IPCW: Cox model for time to surgery on 7 covariates, weight 1/Ŝc(t|X) updated at each censoring time. Surgery-arm IPCW: logistic model for surgery by 182.62 days among those alive, weight 1/p̂ from 182.62 onwards.

**SynPUF (`window_vs_window`, tutorial code 01–02).** For each window w ∈ {30, 90}: clones followed beyond w without having started are censored at w. Weight 1 in [0, w]. After w: logistic model, among clones censored at w and clones that started in the last 7 days of the window, for being uncensored on age, sex and renal disease; starters in the last 7 days get 1/p̂, earlier starters get 1.

**Both.** Weighted Cox model with strategy only, Efron ties, robust (patient-clustered) SE; Wald 95% CI.

## P — Performance measures

Per scenario: bias ± MCSE (including the truth's MC error), empirical SE, model SE and relative error, RMSE, 95% CI coverage, rejection rate of H0: HR = 1 (power in the observed-effect scenario, type I error in the null scenario, and its direction), failures/warnings. **B = 100,000** replicates per scenario for the paper (`b100k` profile; bias MCSE ≈ 0.002 including the truth's MC error from 96 large datasets, coverage MCSE ≈ 0.07%); the default profile (B = 10,000) runs the same design faster. Failures are recorded and excluded; their count is reported.

## Correction (paper Step 7)

- Bias: corrected estimate = observed − bias (the paper's grid search to first order).
- CI: widened by empirical SE / model SE in the observed-effect scenario (when > 1).
- Alternative CI: simulation bootstrap-t, [ψ̂ − q97.5·SE, ψ̂ − q2.5·SE] with q = quantiles of (ψ̂_b − ψ)/SE_b over the observed-effect replicates. Wald, widened and bootstrap-t CIs are compared by coverage and one-sided misses, cross-fitted within each scenario and transferred between scenarios.
- Nominal CI: confidence belt. Extra grid of θ (config/belt.yml): 5 Chebyshev–Lobatto nodes in ψ over the range where CI endpoints fall (observed CI and outer 0.5% of replicate Wald endpoints, widened by 10%), mapped to θ by a monotone spline of ψ(θ) through the calibration and truth values and anchors at θ = ±2 and ±4 (3 large datasets each), limited to the ψ values reached between the anchors; all grid points, including no effect, enter at their CCW large-sample limit ψ(θ) (at no effect ψ(0), not 0), and belt coverage is of ψ(θ); truths from 12 large datasets; replicates per node 1,000, topped up to at most 5,000 until the order-statistic MCSE of both quantiles of T is below 0.05. Quantile curves of T in ψ (weighted cubic polynomial, weights = replicates), CI by test inversion. Checked at θ\* and 0 with that grid point left out (cross-fitted); MCSE of this coverage by a bootstrap over replicates (200 resamples, resampled jointly across grid points), so it includes the Monte Carlo error of the grid quantiles.
- Type I error: p-value = share of null-scenario replicates with |z| ≥ observed |z|; calibrated critical value = 95th percentile of null |z|.
- **Performance of the adjusted analysis (v0.3):** two-fold cross-fitting over replicates (odd/even replicate index): bias, CI inflation and critical value are estimated on one half and applied to the other, in both scenarios; coverage, power and type I error of the adjusted analysis are then averaged over both halves (with two equal halves the adjusted bias is zero by construction, so it is not reported as a result).

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
| 2026-10-24 | Calibration starts from θ̂ and 0 and evaluates each step's large datasets in parallel; evaluations reused for the true values (same seeds). Belt grid: Chebyshev–Lobatto nodes over the needed ψ range instead of a fixed θ list; adaptive replicates per node (target quantile MCSE) | Run time: the belt grid was about two thirds of the computing time and the calibration a serial bottleneck; estimands unchanged |
| 2026-10-25 | Belt: no-effect grid point and coverage check at ψ(0) (CCW limit) instead of 0; nodes placed by a monotone spline with anchors at θ = ±2 and ±4 and limited to the reachable ψ range | Truth 0 at ψ = 0 put a kink in the quantile curves (lopsided misses, observed coverage 94.4%); ψ(θ) levels off for strong effects, so linear extrapolation misplaced the nodes |
| 2026-10-25 | Belt quantile curves: cubic instead of quadratic; the degree no longer reruns the belt simulations | At B = 100,000 the quadratic missed the left-out q_lo at θ* by 0.10 (observed coverage 94.7%, misses 3.1% below / 2.3% above); cubic error ≤ 0.02 |
| 2026-10-26 | Paper reports the b100k run (B = 100,000 per scenario) | Belt coverage check and bias precise enough to be read at the 0.1–0.2% level |
| 2026-10-25 | Belt coverage MCSE by bootstrap including the grid's Monte Carlo error; b100k nodes to quantile MCSE 0.02 (up to 30,000 replicates) | The binomial MCSE (0.07%) ignored the grid quantile error (~0.03 per node), which dominates the check; observed 94.7% / null 95.1% at MCSE 0.035 |
