# model-bias-sim — modelling-bias QBA simulation (CCW case study)

Simulation for Rippin, Luguzis, Antunes, Ralphs: *Extending Quantitative Bias Analysis to Evaluate
Analysis-Induced Bias*
(paper Section 4). Case study: the clone-censor-weighting (CCW) analysis of the Maringe et al. (2020)
lung cancer example (surgery within 6 months vs no surgery).
Adjustment step as in manuscript v0.64 (17 Oct 2026), Section 2.7: inverse calibration (mean CCW estimate
at θ\* = observed estimate), non-parametric bias correction, CI and test, checked in independent simulated
studies, plus a nested (double) simulation that repeats the whole adjustment in each outer dataset.
Design: `protocol/ADEMP.md` (v0.4). Coding standards: `docs/simulation_setup.md`.
How the pipeline fits together (folders, scripts, dependencies): `docs/technical_setup.html`.

## Versions and datasets

| What | Where | Datasets run |
|---|---|---|
| Current analysis for the paper | `main` | Maringe only (`active_datasets` in `config/config.yml`) |
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
Sys.setenv(R_CONFIG_ACTIVE = "smoke"); targets::tar_make()   # smoke run (~30 s)
Sys.unsetenv("R_CONFIG_ACTIVE"); targets::tar_make()         # full run (B = 10,000 + nested)
source("reports/paper/build_section4.R")                     # Word draft of paper Section 4
```

**Run time and parallel runs.** Measured on a 4-core Windows laptop with 3 crew workers (22 Oct 2026):
about 0.17 s per CCW analysis of 200 patients (0.4 s while all cores are busy). The default run takes
about 10–11 hours elapsed, roughly 30 CPU-hours: inverse calibration 1 hour (sequential, 4 steps of 2,000
studies), nested simulation about 27 CPU-hours (outer datasets at θ\* ~130 s each, with no effect ~65 s
each), the 2 × 10,000 main and 2 × 10,000 check studies about 2 CPU-hours, true values well under 1.
With `crew` installed, branches run on parallel local workers (default 2, `workers` in
`config/config.yml`); `Sys.setenv(QBA_WORKERS = parallel::detectCores() - 1)` before `tar_make()` uses
more cores without editing the config. Results do not depend on the number of workers.
Do not edit `config/config.yml` during or between runs unless you want to rerun everything: every target
depends on it. For the CCW and simple-adjustment results without waiting for the nested simulation:
`targets::tar_make(names = c("adjustment_maringe", "validation_maringe"))`; a later full `tar_make()`
keeps everything already built, including finished nested branches.

**More replicates.** `B` (replicates per scenario) is in `config/config.yml`. For a precision check
without touching the paper run, use the `b100k` profile (100,000 main and 100,000 check studies per
scenario, 96 instead of 24 large datasets for the true value, nested 1,000 outer × 500 inner) with its own
store;
figures, report and Word draft go to `reports/b100k/` (the paper outputs in `reports/` are untouched):

```r
Sys.setenv(R_CONFIG_ACTIVE = "b100k", QBA_WORKERS = parallel::detectCores() - 1)
targets::tar_make(store = "_targets_b100k")      # store name must be _targets_<profile>
source("reports/paper/build_section4.R")         # -> reports/b100k/paper/section4_draft.docx
targets::tar_read(performance_maringe, store = "_targets_b100k")
Sys.unsetenv("R_CONFIG_ACTIVE")                  # back to the paper run
```

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
| `R/` | data, dgm, truth (incl. inverse calibration), pipeline, simulate, performance, correction (simple adjustment and its check), nested, plots, paper, belt (optional), utils |
| `config/config.yml` | seeds, B, check and nested sizes, truth size, calibration, datasets (`default`, `smoke`, `two_datasets`, `b100k`) |
| `config/belt.yml` | optional bootstrap-t / confidence-belt comparison (`enabled: false` by default) |
| `data/raw/` | Maringe lung cancer CSV and SynPUF tutorial cohort (read-only) + source/checksums |
| `tests/testthat/` | unit tests |
| `reports/` | paper outputs, protocol figure, Quarto report (two-dataset profile) |
| `protocol/ADEMP.md` | simulation protocol and change log |
| `docs/` | coding standards and technical overview |
| `archive/` | frozen results of earlier runs |
