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

#' Software used for the replicates (reported with the run time; the number of workers is
#' counted from the timestamps by `observed_workers()`)
#' @param sims Replicate results; passed so this is re-recorded only when the replicates rerun.
#' @return List: crew version, cores available, number of replicate rows.
record_run_settings <- function(sims) {
  list(
    crew = if (requireNamespace("crew", quietly = TRUE)) {
      as.character(utils::packageVersion("crew"))
    } else {
      "not installed"
    },
    cores = parallel::detectCores(),
    replicates = nrow(sims)
  )
}

#' Number of parallel workers a run actually used, read from the targets metadata
#' @param meta `tar_meta()` rows with name, time (finish) and seconds (run time).
#' @return Integer: highest number of simulation branches running at the same time for a
#'   sustained period (at least 30 s or 5% of the stage); brief overlaps from the slight delay
#'   in recording finish times are ignored.
observed_workers <- function(meta) {
  # Branches of the simulation patterns (branch-name hashes differ between targets versions)
  is_branch <- if ("type" %in% names(meta)) meta$type == "branch" else TRUE
  m <- meta[is_branch & grepl("^(sim_results|truth_rows|calib_rows_[0-9]|belt_truth_rows|belt_sims_[0-9])_",
                              meta$name), ]
  if (nrow(m) == 0) return(NA_integer_)
  st <- as.numeric(m$time) - m$seconds
  en <- as.numeric(m$time)
  ev <- data.frame(t = c(st, en), d = c(rep(1L, length(st)), rep(-1L, length(en))))
  ev <- ev[order(ev$t, ev$d), ]
  level <- cumsum(ev$d)
  held <- tapply(c(diff(ev$t), 0), level, sum)
  span <- max(en) - min(st)
  lv <- as.integer(names(held))[held >= max(30, 0.05 * span)]
  max(lv[lv > 0], 1L)
}
