# qba-sim — modelling-bias QBA simulation (CCW case study)

Simulation for Rippin, Luguzis, Antunes, Ralphs: *Quantifying modelling bias in pharmacoepidemiology*.
Design: `protocol/ADEMP.md`. Coding standards: `simulation_setup.md` (project doc).

## Reproduce

```r
# once, on your machine (records your package versions in renv.lock)
install.packages("renv")
renv::init(bare = TRUE)
renv::install(c("targets", "tarchetypes", "crew", "config", "data.table", "survival", "haven",
                "checkmate", "cli", "withr", "rsimsum", "ggplot2", "scales", "quarto",
                "sessioninfo", "testthat"))
renv::snapshot()

# every run
testthat::test_dir("tests/testthat")            # unit tests (~10 s)
Sys.setenv(R_CONFIG_ACTIVE = "smoke"); targets::tar_make()   # smoke run (~30 s)
Sys.unsetenv("R_CONFIG_ACTIVE"); targets::tar_make()         # full run (~20 min, 2 workers)
```

The report is `reports/qba_report.html` (rendered by the pipeline from the targets store).
Results of the 5 Oct 2026 run (two datasets): unzip `_targets_store.zip` here to get `_targets/`; then `targets::tar_read(performance)` etc. work without rerunning.
Quarto CLI must be installed for the report target.

## Layout

| Path | Content |
|---|---|
| `_targets.R` | pipeline (only entry point) |
| `R/` | data, dgm, truth, pipeline, simulate, performance, correction, plots, utils |
| `config/config.yml` | seeds, B, grid, truth size, pipeline settings (`default`, `smoke`) |
| `data/raw/` | Maringe lung cancer CSV and SynPUF tutorial cohort (read-only) + source/checksums |
| `tests/testthat/` | unit tests |
| `reports/` | Quarto report and figures |

Runs of 2 and 5 Oct 2026 used R 4.3.3, targets 1.4.1, crew 0.8.0 (Linux, 2 workers).
