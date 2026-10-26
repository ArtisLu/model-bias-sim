# model-bias-sim — modelling-bias QBA simulation (CCW case study)

Simulation for Rippin, Luguzis, Antunes, Ralphs: *Quantifying modelling bias in pharmacoepidemiology*
(paper Section 4). Case study: the clone-censor-weighting (CCW) analysis of the Maringe et al. (2020)
lung cancer example (surgery within 6 months vs no surgery).
Design: `protocol/ADEMP.md`. Coding standards: `docs/simulation_setup.md`.
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
Sys.setenv(R_CONFIG_ACTIVE = "b100k")                        # paper run (B = 100,000)
targets::tar_make(store = "_targets_b100k")
source("reports/paper/build_section4.R")                     # Word draft of paper Section 4
```

**Parallel runs.** With `crew` installed, the replicates and the true values run on parallel local
workers (default 2, `workers` in `config/config.yml`). To use more cores without editing the config,
set `Sys.setenv(QBA_WORKERS = parallel::detectCores() - 1)` before `tar_make()`. Results do not depend
on the number of workers. Each large dataset behind a true value is its own branch, so calibration
steps, true values and confidence-belt truths also spread over the workers; only the 2–4 calibration
steps follow each other. Large datasets need ~1.5 GB each, so the number of workers is limited by memory.

**Paper run and profiles.** The paper reports the `b100k` profile: 100,000 replicates per scenario,
96 large datasets per true value, and confidence-belt grid points topped up to a quantile MCSE of 0.02
(up to 30,000 replicates each). It has its own store; figures, report and Word draft go to
`reports/b100k/`. The default profile (B = 10,000, outputs in `reports/`) is a faster run of the same
design. `B` and the other settings are in `config/config.yml` and `config/belt.yml`.

```r
Sys.setenv(R_CONFIG_ACTIVE = "b100k", QBA_WORKERS = 3)  # ~1.5 GB RAM per worker
targets::tar_make(store = "_targets_b100k")      # store name must be _targets_<profile>
source("reports/paper/build_section4.R")         # -> reports/b100k/paper/section4_draft.docx
targets::tar_read(performance_maringe, store = "_targets_b100k")
Sys.unsetenv("R_CONFIG_ACTIVE")                  # back to the default profile
```

The project folder may be synced (Dropbox): the report is rendered in a temporary folder for that
reason. Excluding the `_targets*` stores from syncing avoids file locks and large uploads.

Results of the 12 Oct 2026 run (Maringe, B = 10,000, R 4.3.3, 2 crew workers, UTF-8 locale): unzip
`results/targets_store_2026-10-12_maringe.zip` here to get `_targets/`. `targets::tar_outdated()` is then
empty and `targets::tar_read()` works without rerunning (a full run takes about 21 minutes with 2 workers,
34 minutes sequentially). The store does not depend on the project folder, but it does depend on the
locale: run R in a UTF-8 locale (default on Windows with R >= 4.2).

Outputs: `reports/paper/<dataset>/` (paper figures as PNG and EPS), `reports/paper/section4_draft.docx`
and `reports/fig_protocol_<dataset>.png`. The report `reports/qba_report.html` is rendered by the pipeline
from `reports/qba_report.qmd` when Maringe is active and the Quarto CLI is installed (not in smoke runs).

## Layout

| Path | Content |
|---|---|
| `_targets.R` | pipeline (only entry point) |
| `R/` | data, dgm, truth, pipeline, simulate, performance, correction, plots, paper, utils |
| `config/config.yml` | seeds, B, truth size, datasets, pipeline and censoring settings (`default`, `smoke`, `two_datasets`) |
| `data/raw/` | Maringe lung cancer CSV and SynPUF tutorial cohort (read-only) + source/checksums |
| `tests/testthat/` | unit tests |
| `reports/` | paper outputs, protocol figure, Quarto report (two-dataset profile) |
| `protocol/ADEMP.md` | simulation protocol and change log |
| `docs/` | coding standards and technical overview |
| `archive/` | frozen results of earlier runs |
