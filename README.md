# model-bias-sim — modelling-bias QBA simulation (CCW case study)

Simulation for Rippin, Luguzis, Antunes, Ralphs: *Extending Quantitative Bias Analysis to Evaluate
Analysis-Induced Bias*
(paper Section 4). Case study: the clone-censor-weighting (CCW) analysis of the Maringe et al. (2020)
lung cancer example (surgery within 6 months vs no surgery).
Adjustment step as in manuscript v0.65 (24 Oct 2026), Section 2.7: inverse calibration (mean CCW estimate
at θ\* = observed estimate, search started at θ = −0.7), non-parametric bias correction, CI and test,
checked in independent simulated studies, plus a nested (double) simulation that repeats the whole
adjustment in each outer dataset. Additional analyses of Section 4.6 (second endpoint: difference in
12-month survival; generalized-gamma prediction model; hypothetical-trial target for an estimator that
is not consistent) and the pilot study for the outer / inner split (Appendix). Performance measures:
bias, empirical and model SE, RMSE, MSE (= RMSE²), coverage, type I error, each with MCSE; no power.
Design: `protocol/ADEMP.md` (v0.5). Coding standards: `docs/simulation_setup.md`.
How the pipeline fits together (folders, scripts, dependencies): `docs/technical_setup.html`.

## Versions and datasets

| What | Where | Datasets run |
|---|---|---|
| Current analysis for the paper | `main` | Maringe: log HR (`maringe`) and 12-month survival difference (`maringe_surv12`) (`active_datasets` in `config/config.yml`) |
| Two-dataset comparison of 5 Oct 2026 (frozen) | branch `archive/two-datasets` | Maringe + SynPUF |
| Results of that run | `archive/2026-10-05_two_datasets/` | report, figures, targets store |

The SynPUF code is still on `main` and covered by the unit tests and the smoke run; it is just not run
by default.
To rerun the two-dataset comparison with the current code:
`Sys.setenv(R_CONFIG_ACTIVE = "two_datasets"); targets::tar_make()`.

## Reproduce

```r
# once, on your machine (records your package versions in renv.lock)
install.packages("renv")
renv::init(bare = TRUE)
renv::install(c("targets", "tarchetypes", "crew", "config", "data.table", "survival", "haven",
                "checkmate", "cli", "withr", "rsimsum", "ggplot2", "scales", "quarto",
                "sessioninfo", "testthat", "officer", "flextable"))
renv::snapshot()

# every run
testthat::test_dir("tests/testthat")                         # unit tests (~10 s)
Sys.setenv(R_CONFIG_ACTIVE = "smoke"); targets::tar_make()   # smoke run (~1 min)
Sys.setenv(R_CONFIG_ACTIVE = "b100k")                        # paper run (B = 100,000 + nested)
targets::tar_make(store = "_targets_b100k")
source("reports/paper/build_section4.R")                     # Word draft of paper Section 4
```

**Run time and parallel runs.** Measured on a 4-core Windows laptop with 3 crew workers (22 Oct 2026):
about 0.17 s per CCW analysis of 200 patients (0.4 s while all cores are busy). The default run takes
about 10–11 hours elapsed, roughly 30 CPU-hours: inverse calibration 1 hour (sequential, 4 steps of 2,000
studies), nested simulation about 27 CPU-hours (outer datasets at θ\* ~130 s each, with no effect ~65 s
each), the 2 × 10,000 main and 2 × 10,000 check studies about 2 CPU-hours, true values well under 1.
With `crew` installed, branches run on parallel local workers (default 2, `workers` in
`config/config.yml`); `Sys.setenv(QBA_WORKERS = parallel::detectCores() - 1)` before `tar_make()` uses
more cores without editing the config. Results do not depend on the number of workers. Each large
dataset behind a true value (and behind the optional belt) is its own branch, so truths spread over the
workers; large datasets need ~1.5 GB each, so the number of workers is limited by memory.
Do not edit `config/config.yml` during or between runs unless you want to rerun everything: every target
depends on it. For the CCW and simple-adjustment results without waiting for the nested simulation:
`targets::tar_make(names = c("adjustment_maringe", "validation_maringe"))`; a later full `tar_make()`
keeps everything already built, including finished nested branches.

**Pilot study and outer / inner split (paper Appendix).** The targets `pilot_nested_<ds>`, `pilot_timing_<ds>`,
`allocation_<ds>` and `pilot_anova_<ds>` run `pilot:` (50 outer × 50 inner by default) and give k = σ_w²/V,
the costs of an outer replicate and of an inner run, R\*, the outer runs for the target precision and the
cost of the configured design (`R/allocation.R`). Run them first
(`targets::tar_make(names = starts_with("allocation_"))`), then set `nested: outer` and `inner` in
`config/config.yml` for the paper run. The configured values are the result of the pilot of 27 Oct 2026
(see the change log in `protocol/ADEMP.md`). `nested: shrink: true` applies the shrinkage of the outer
deviations by √(V/Var_obs).

**Additional analyses (paper 4.6).**

```r
# Generalized-gamma prediction model (own store; the paper run is untouched)
Sys.setenv(R_CONFIG_ACTIVE = "gengamma", QBA_WORKERS = parallel::detectCores() - 1)
targets::tar_make(store = "_targets_gengamma")
Sys.unsetenv("R_CONFIG_ACTIVE")
```

The 12-month survival endpoint is a second dataset entry (`maringe_surv12`) with the same data and prediction
model. The hypothetical-trial target (`truth_rct_*`, `adjustment_rct_*`, `nested_rct_*`, `perf_compare_rct_*`;
HR endpoint only; assumptions in `rct:` and `protocol/ADEMP.md`) reuses the replicates and the inner runs;
only the truth changes.

**Running on a high-performance server.** Nothing in the code is tied to a machine: `QBA_WORKERS` sets the
number of crew workers (default `workers` in the config), branches are independent (common random numbers
tied to the replicate index), and results do not depend on the number of workers. To move the project:
copy the repository (without `_targets*` stores), `renv::restore()` from `renv.lock` (R 4.3.3 or newer;
install the system libraries of the packages first), run the unit tests and the smoke profile, then
`Sys.setenv(QBA_WORKERS = <cores - 1>); targets::tar_make()` in a persistent session (`tmux`, `screen` or the
scheduler of the server; for a SLURM cluster `crew.cluster::crew_controller_slurm()` can replace the local
controller in `_targets.R`). Memory: about 1.5 GB per worker for the large-sample truth (50,000 patients
for the main truth, 25,000 in the nested runs); the stores are written to the project folder, so use fast
local disk.

**Paper run and profiles.** The paper reports the `b100k` profile: 100,000 replicates per scenario,
96 large datasets per true value, and the nested simulation with 2,000 outer × 50 inner runs (chosen by the
pilot, see below). It has its own store; figures, report and Word draft go to `reports/b100k/`. The
default profile (B = 10,000, outputs in `reports/`) is a faster run of the same design. `B` and the other
settings are in `config/config.yml`. The confidence belt of the earlier design (adaptive grid, config in
`config/belt.yml`, cubic quantile curves) is optional since manuscript v0.64: `enabled: true` in
`config/belt.yml` builds it next to the non-parametric adjustment (the smoke profile does).

```r
Sys.setenv(R_CONFIG_ACTIVE = "b100k", QBA_WORKERS = 3)  # ~1.5 GB RAM per worker
targets::tar_make(store = "_targets_b100k")      # store name must be _targets_<profile>
source("reports/paper/build_section4.R")         # -> reports/b100k/paper/section4_draft.docx
targets::tar_read(performance_maringe, store = "_targets_b100k")
Sys.unsetenv("R_CONFIG_ACTIVE")                  # back to the default profile
```

The project folder may be synced (Dropbox): the report is rendered in a temporary folder for that
reason. Excluding the `_targets*` stores from syncing avoids file locks and large uploads.

Results of the 12 Oct 2026 run (protocol v0.3: calibration on the true value, cross-fitted adjustment;
Maringe, B = 10,000, R 4.3.3, 2 crew workers, UTF-8 locale): unzip
`results/targets_store_2026-10-12_maringe.zip` into a separate checkout of commit `31508f6` to read it with
`targets::tar_read()`; the v0.4 code has different targets. The store does not depend on the project folder, but it does depend on the
locale: run R in a UTF-8 locale (default on Windows with R >= 4.2).

Outputs: `reports/paper/<dataset>/` (paper figures as PNG and EPS: performance, distributions, nested
simulation, model fit, balance), `reports/paper/section4_draft.docx`
and `reports/fig_protocol_<dataset>.png`. The report `reports/qba_report.html` is rendered by the pipeline
from `reports/qba_report.qmd` when Maringe is active and the Quarto CLI is installed (not in smoke runs).

## Layout

| Path | Content |
|---|---|
| `_targets.R` | pipeline (only entry point) |
| `R/` | data, dgm (Weibull, generalized gamma), truth (incl. inverse calibration), pipeline (log HR, 12-month survival difference), simulate, performance, correction (simple adjustment and its check), nested, allocation (pilot, outer / inner split), rct (hypothetical-trial target), plots, paper, belt (optional), utils |
| `config/config.yml` | seeds, B, check, nested and pilot sizes, truth size, calibration, `rct`, datasets (profiles `default`, `smoke`, `two_datasets`, `gengamma`, `b100k`) |
| `config/belt.yml` | optional bootstrap-t / confidence-belt comparison (`enabled: false` by default) |
| `data/raw/` | Maringe lung cancer CSV and SynPUF tutorial cohort (read-only) + source/checksums |
| `tests/testthat/` | unit tests |
| `reports/` | paper outputs, protocol figure, Quarto report (two-dataset profile) |
| `protocol/ADEMP.md` | simulation protocol and change log |
| `docs/` | coding standards and technical overview |
| `archive/` | frozen results of earlier runs |
