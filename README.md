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
Sys.unsetenv("R_CONFIG_ACTIVE"); targets::tar_make()         # full run (B = 10,000)
source("reports/paper/build_section4.R")                     # Word draft of paper Section 4
```

Results of the 12 Oct 2026 run (Maringe, B = 10,000, R 4.3.3, UTF-8 locale): unzip
`results/targets_store_2026-10-12_maringe.zip` here to get `_targets/`; `targets::tar_read()` then works
without rerunning (about 35 minutes on one core).

Outputs: `reports/paper/<dataset>/` (paper figures as PNG and EPS), `reports/paper/section4_draft.docx`
and `reports/fig_protocol_<dataset>.png`. The Quarto colleague report (`reports/qba_report.qmd`) is
written for two datasets and runs only in the `two_datasets` profile.

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
