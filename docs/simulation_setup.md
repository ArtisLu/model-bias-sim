# Simulation set-up and coding standards — Modelling-bias QBA

**Project:** Quantifying modelling bias in pharmacoepidemiology (Rippin, Ralphs, Luguzis et al.)
**Owner of the simulation:** Artis Luguzis
**Version:** 1.0 (2 Oct 2026)
**Purpose of this document:** the standing reference for *how* the simulation is built, run and reported. The *content* of the simulation (aims, DGM, estimand, pipeline, grid) lives in the ADEMP protocol (`protocol/ADEMP.md`). Any R code written for this project must follow this document unless the protocol says otherwise.

> **Current scope (Sep 2026):** one base dataset (a public CCW dataset), one analysis pipeline (CCW), one scenario. Truth = the treatment parameter of the CCW model. The set-up below is generic so that further pipelines or scenarios can be added as rows in the scenario table, without changing the code structure.

---

## 1. Guiding principles

1. **Protocol first.** The design is written down (ADEMP, Morris, White & Crowther, *Stat Med* 2019) before any production run. Changes after the first full run are logged in the protocol's change log, with reason and date.
2. **Reproducible.** Anyone with the repository runs one command and gets identical numbers: pinned package versions, explicit seeds, no manual steps, no hidden state.
3. **Replicable.** The protocol and code are enough for an independent team to re-implement the study and reach the same conclusions. The code is archived with a DOI at submission.
4. **The pipeline under test is one function.** The same function object estimates the effect on the real data and on every simulated dataset. Nothing is copied or re-typed.
5. **Separate the stages.** Data → DGM → truth → simulate → pipeline → performance → correction. Each stage is its own function and its own pipeline target, with outputs cached.
6. **Store raw replicate results.** Every replicate's estimate, SE, CI, fit status and seed is saved. Performance measures are always derived from the stored results, never computed on the fly and discarded.
7. **Quantify Monte Carlo error.** Every performance measure is reported with its Monte Carlo SE (MCSE). The number of replicates B is justified by a target MCSE.
8. **Failures are data.** Non-converged or failed fits are recorded, counted and reported, never silently dropped. The rule for handling them is pre-specified.
9. **Verify small before running big.** Unit tests and a smoke run (small B, few grid points) pass before every full run.

---

## 2. Software stack

| Purpose | Tool | Notes |
|---|---|---|
| R version | Current R release | Recorded in `renv.lock`; do not upgrade mid-study |
| Package versions | **renv** | `renv::snapshot()` after any package change; lockfile is committed |
| Pipeline / caching | **targets** + **tarchetypes** | `tar_make()` is the only entry point; reruns only what changed |
| Parallel execution | **crew** (`crew_controller_local`) | Parallelises across targets/branches; results must not depend on worker count |
| Data handling | **data.table** | For replicate-level tables; native pipe `\|>` elsewhere |
| Survival / CCW | **survival**, **flexsurv** (or **rstpm2**), **WeightIt**; CCW per chosen dataset tutorial (e.g. **ccwr**) | Pipeline packages are fixed by the protocol |
| Missing data (if used) | **mice** | Seed passed explicitly |
| Performance measures | **rsimsum** | Bias, empSE, modSE, coverage, power, all with MCSE |
| Seeds | **withr** (`with_seed()`), L'Ecuyer-CMRG RNG | See Section 5 |
| Configuration | **config** (YAML) | All tunables in `config/config.yml`; no magic numbers in functions |
| Argument checks | **checkmate** | Assert inputs at the top of each exported function |
| Messages | **cli** | No `print()`/`cat()` for progress |
| Tests | **testthat** | `tests/testthat/`; run by `testthat::test_dir("tests/testthat")` |
| Formatting / linting | **Air** (formatter), **lintr** | `air.toml` and `.lintr` in repo; tidyverse style |
| Figures / tables | **ggplot2**, **gt** | Saved by targets, not by hand |
| Report | **Quarto** (`tarchetypes::tar_quarto()`) | Report reads results with `tar_read()`; contains no computation |
| Session record | **sessioninfo** | Stored as a target with each run |
| Version control | **git** (+ GitHub/IQVIA repo) | Commit before every full run; archive release on Zenodo at submission |

Optional for full computational reproducibility: a **rocker** Docker image with the same R version (useful if runs move to a server).

---

## 3. Project structure

```
qba-sim/
├── _targets.R            # pipeline definition (the only entry point)
├── R/                    # functions only, no top-level code with side effects
│   ├── data.R            # load and prepare the base dataset
│   ├── dgm.R             # fit reference DGM; simulate outcomes
│   ├── truth.R           # compute the true value psi(theta)
│   ├── pipeline.R        # the analysis pipeline under test
│   ├── simulate.R        # replicate/batch runner, safe wrapper
│   ├── performance.R     # performance measures with MCSE
│   ├── correction.R      # bias curve, inversion, CI and type I calibration
│   ├── plots.R           # figure functions
│   └── utils.R           # seeds, helpers
├── config/config.yml     # seeds, B, grid, batch size, workers, pipeline settings
├── protocol/ADEMP.md     # design + change log
├── data/raw/             # base data as obtained (read-only; never edited)
├── tests/testthat/       # unit tests
├── reports/qba_report.qmd
├── renv.lock, .Rprofile, renv/
├── air.toml, .lintr, .gitignore   # _targets/ store is git-ignored
└── README.md             # how to reproduce: renv::restore(); targets::tar_make()
```

---

## 4. Pipeline (`_targets.R`)

Skeleton only; target names are fixed, function bodies follow the protocol.

```r
# _targets.R — modelling-bias QBA simulation. Run with targets::tar_make().
library(targets)
library(tarchetypes)

tar_source("R")
cfg <- config::get(file = "config/config.yml")

tar_option_set(
  packages   = c("data.table", "survival", "checkmate"),
  seed       = cfg$seed_master,
  controller = crew::crew_controller_local(workers = cfg$workers)
)

list(
  # Data and observed analysis
  tar_target(base_file, cfg$base_data, format = "file"),
  tar_target(base_data, prepare_base_data(base_file)),
  tar_target(obs_result, run_pipeline(base_data, cfg$pipeline)),

  # Simulated world
  tar_target(dgm, fit_dgm(base_data, cfg$dgm)),
  tar_target(theta_grid, make_grid(obs_result, cfg$grid)),
  tar_target(truth, compute_truth(dgm, base_data, theta_grid), pattern = map(theta_grid)),
  tar_target(crn_seeds, make_crn_seeds(cfg$seed_master, cfg$B)),
  tar_target(batch_id, seq_len(cfg$B %/% cfg$batch_size)),

  # Replicates: one branch per (theta, batch)
  tar_target(
    sim_results,
    run_batch(base_data, dgm, theta_grid, crn_seeds, batch_id, cfg),
    pattern = cross(theta_grid, batch_id)
  ),

  # Summaries and correction
  tar_target(performance, summarise_performance(sim_results, truth)),
  tar_target(bias_curve, fit_bias_curve(performance)),
  tar_target(correction, correct_estimate(obs_result, bias_curve, performance)),
  tar_target(verification, verify_correction(base_data, dgm, correction, crn_seeds, cfg)),

  # Outputs
  tar_target(fig_bias_curve, plot_bias_curve(performance, bias_curve)),
  tar_target(session, sessioninfo::session_info()),
  tar_quarto(report, "reports/qba_report.qmd")
)
```

Rules:

- Replicates are **batched** (e.g. 50–100 per branch) to keep target overhead low; one branch = one (θ, batch) pair.
- Store **summaries only** (one row per replicate), never fitted model objects.
- A smoke run uses a separate config profile (`R_CONFIG_ACTIVE=smoke`: small B, 2–3 grid points, 1–2 workers).
- Interrupted runs resume automatically; do not delete `_targets/` to "start clean" unless the protocol changed.

---

## 5. Random numbers and seeds

**Requirements**

1. One master seed in `config.yml`. All other seeds are derived from it deterministically.
2. **Common random numbers (CRN):** replicate *b* uses the same draws at every grid value θ. Therefore draws are tied to the replicate index *b*, not to the target, branch or worker. (Do not rely on targets' automatic per-target seeds for outcome draws: they differ between branches, which breaks CRN.)
3. Results are identical for any number of workers and any batch size. This is tested (Section 8).
4. Event times are drawn by **inversion** of the DGM's cumulative hazard, T = H⁻¹(−log U | X, θ), so that the same U gives monotone, smoothly varying times across θ.
5. Each random component (event times, censoring, MI, bootstrap) draws from the replicate's seed in a **fixed order**. Adding a component appends to that order; it never reorders it.
6. The replicate seed is stored with each result row.

**Template**

```r
# Seeds for replicates 1..B, derived once from the master seed.
make_crn_seeds <- function(seed_master, B) {
  withr::with_seed(
    seed_master,
    sample.int(.Machine$integer.max, B),
    .rng_kind = "L'Ecuyer-CMRG"
  )
}

# All random draws for one replicate, in a fixed order.
draw_replicate_rng <- function(seed, n) {
  withr::with_seed(seed, .rng_kind = "L'Ecuyer-CMRG", {
    list(
      u_event  = runif(n),
      u_censor = runif(n),
      seed_aux = sample.int(.Machine$integer.max, 1) # for MI / bootstrap
    )
  })
}
```

---

## 6. Function design and coding style

**Structure**

- Every stage is a **pure function**: inputs as arguments, output as return value. No global variables, no `rm(list = ls())`, no `setwd()`, no `set.seed()` outside the seed helpers.
- Verb-based `snake_case` names: `fit_dgm()`, `simulate_dataset()`, `compute_truth()`, `run_pipeline()`.
- Settings come in through a `spec`/`cfg` list from `config.yml`; defaults live in the config, not in function signatures.
- Validate inputs with `checkmate` at the top of each function.
- Functions do one thing and fit on one screen; split otherwise.
- Use `pkg::fun()` for packages not in `tar_option_set(packages = ...)`.
- Within a function, use one data paradigm (data.table *or* base/tidyverse), not a mix.

**Pipeline function contract**

`run_pipeline(data, spec)` returns a one-row `data.table` with fixed columns:
`estimate, se, ci_lo, ci_hi, p_value, n_events, status, message`.
The same function runs on the real data and on every simulated dataset.

**Failure handling**

Warnings are captured, not treated as failures by default (dropping fits on warning is itself a source of bias). What counts as a failure, and whether failed replicates are excluded or retried, is pre-specified in the protocol.

```r
# Run the pipeline, recording errors and warnings instead of stopping.
safely_run <- function(fun, ...) {
  warns <- character()
  res <- withCallingHandlers(
    tryCatch(fun(...), error = function(e) e),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (inherits(res, "error")) {
    return(data.table::data.table(
      estimate = NA_real_, se = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
      p_value = NA_real_, n_events = NA_integer_,
      status = "error", message = conditionMessage(res)
    ))
  }
  if (length(warns)) {
    res[, `:=`(status = "warning", message = paste(unique(warns), collapse = " | "))]
  }
  res
}
```

**Comments and documentation**

- Each file starts with a 1–3 line header: what the file contains.
- Each function has a short roxygen-style header: one-line description, key `@param`s, `@return`. No tutorials.
- Inline comments explain *why* (a modelling choice, a protocol reference such as `# ADEMP D2`), not *what* the code does.
- Formatting is done by **Air** (run before every commit); `lintr` must show no warnings.

---

## 7. Monte Carlo precision and performance

- **Choose B** from a pilot: B ≈ (empirical SD / target MCSE)². For coverage, MCSE = √(p(1 − p)/B); a 0.5% MCSE at p = 0.95 needs ≈ 1,900 replicates. Record the targets and the resulting B in the protocol.
- **Measures per grid point:** bias ± MCSE, standardised bias (bias / empSE), empirical SE, model SE and its relative error, RMSE, coverage, rejection rate (type I error at θ = 0, power elsewhere), failure rate.
- Use `rsimsum::simsum()` on the stored replicate table; truth is joined by θ from the `truth` target.
- **Relevance thresholds** for bias (on the log scale) are pre-specified. With large B, p-values for bias are not used.
- Report bias on the log scale; exp(bias) − 1 is the relative bias on the ratio scale. Never as % of log-HR.

---

## 8. Verification and testing

Run before every full run (`testthat::test_dir("tests/testthat")` + smoke profile):

| Test | Pass criterion |
|---|---|
| DGM recovers its own parameters on a large simulated dataset | Within tolerance set in test |
| Truth: analytic vs large-sample Monte Carlo (where both exist) | Agree within tolerance |
| Pipeline on real data equals `obs_result` | Identical |
| Same seed → identical replicate | Identical |
| CRN: replicate *b* uses the same `u_event` at every θ | Identical |
| Results independent of workers (1 vs 4) and batch size | Identical |
| At θ such that the DGM equals the analysis model (if applicable) | No bias beyond MCSE |
| `safely_run()` returns a valid row for an error and a warning | Correct status |

---

## 9. Reproducibility checklist (each full run)

- [ ] Working tree committed; git hash stored in results metadata.
- [ ] `renv::status()` clean; `renv.lock` committed.
- [ ] Tests and smoke run pass.
- [ ] `tar_make()` from a clean session; no manual steps.
- [ ] `session` target updated; report rendered from the targets store.
- [ ] Protocol change log updated if anything in the design changed.
- [ ] Base data: source, version and checksum recorded (`data/raw/README`).

---

## 10. Reporting

- Follow the ADEMP reporting items of Morris et al. (2019): aims, DGM, estimands, methods, performance measures, B and its justification, MCSE, software versions, failure handling.
- Main text: one worked end-to-end example (observed → bias curve → corrected estimate and CI → verification).
- Supplement: full performance tables, failure rates, sensitivity to the DGM, code and session information.
- All tables and figures are produced by the pipeline; none are edited by hand.

---

## Reference instruction (paste into chat or project instructions)

> For any R code in this project, follow `simulation_setup.md`: targets + crew pipeline, renv, config in `config/config.yml`, modular pure functions in `R/` by stage (data → DGM → truth → simulate → pipeline → performance → correction), CRN seeds tied to replicate index, per-replicate results stored, failures recorded, MCSE reported, Air formatting, short roxygen headers and brief why-comments only.
