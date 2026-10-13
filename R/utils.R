# Seeds, common random numbers and the safe pipeline wrapper.

#' Seeds for replicates 1..B, derived once from the master seed
#' @param seed_master Integer master seed from config.
#' @param B Number of replicates.
#' @return Integer vector of length B.
make_crn_seeds <- function(seed_master, B) {
  checkmate::assert_count(seed_master)
  checkmate::assert_count(B, positive = TRUE)
  withr::with_seed(
    seed_master,
    sample.int(.Machine$integer.max, B),
    .rng_kind = "L'Ecuyer-CMRG"
  )
}

#' All random draws for one replicate, in a fixed order
#' @param seed Replicate seed.
#' @param n Number of patients.
#' @return List with `u_event` (event-time inversion), `seed_aux`, `u_cens` (censoring times).
draw_replicate_rng <- function(seed, n) {
  withr::with_seed(seed, .rng_kind = "L'Ecuyer-CMRG", {
    list(
      u_event = runif(n),
      # Reserved for later components (MI, bootstrap); appended, never reordered
      seed_aux = sample.int(.Machine$integer.max, 1),
      u_cens = runif(n)
    )
  })
}

#' Empty result row in the pipeline contract format
#' @param status,message Failure status and message.
#' @return One-row data.table.
empty_result <- function(status, message) {
  data.table::data.table(
    estimate = NA_real_, se = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
    p_value = NA_real_, n_events = NA_integer_,
    status = status, message = message
  )
}

#' Run the pipeline, recording errors and warnings instead of stopping
#' @param fun Pipeline function.
#' @param ... Arguments passed to `fun`.
#' @return One-row data.table in the pipeline contract format.
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
    return(empty_result("error", conditionMessage(res)))
  }
  if (length(warns)) {
    res[, `:=`(status = "warning", message = paste(unique(warns), collapse = " | "))]
  }
  res
}

#' How the replicates were run (reported with the run time)
#' @param default_workers `cfg$workers`; QBA_WORKERS overrides it as in `_targets.R`.
#' @param sims Replicate results; passed so this is re-recorded only when the replicates rerun.
#' @return List: workers, crew version, cores, number of replicate rows.
record_run_settings <- function(default_workers, sims) {
  # Read at run time, not from a global, so rebuilding only the report cannot change it
  w <- as.integer(Sys.getenv("QBA_WORKERS", default_workers))
  crew_ok <- w > 1 && requireNamespace("crew", quietly = TRUE)
  list(
    workers = if (crew_ok) w else 1L,
    crew = if (crew_ok) as.character(utils::packageVersion("crew")) else "not installed",
    cores = parallel::detectCores(),
    replicates = nrow(sims)
  )
}
