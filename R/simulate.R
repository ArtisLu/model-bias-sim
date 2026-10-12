# Replicate/batch runner: simulate outcomes and run the pipeline (paper Steps 3-5).

#' Run one batch of replicates at one theta
#' @param data Base data; @param dgm Fitted DGM; @param theta Scenario value (scalar).
#' @param seeds CRN seeds (length B); @param batch_id Batch index.
#' @param pipeline Pipeline spec; @param batch_size Replicates per batch.
#' @return data.table, one row per replicate (pipeline contract + theta, rep, seed).
run_batch <- function(data, dgm, theta, seeds, batch_id, pipeline, batch_size) {
  checkmate::assert_number(theta)
  checkmate::assert_count(batch_id, positive = TRUE)
  reps <- (batch_id - 1L) * batch_size + seq_len(batch_size)
  reps <- reps[reps <= length(seeds)]
  out <- lapply(reps, function(b) {
    # Draws depend on the replicate index only, so replicate b is identical in both scenarios
    rng <- draw_replicate_rng(seeds[b], nrow(data))
    sim <- simulate_outcomes(data, dgm, theta, rng$u_event, rng$u_cens)
    res <- safely_run(run_pipeline, sim, pipeline)
    res[, `:=`(
      theta = theta, rep = b, seed = seeds[b],
      n_outcomes = sum(sim$status), n_treated = sum(is.finite(sim$trt_time)),
      n_lost = sum(sim$status == 0L & sim$time < sim$admin_time)
    )]
  })
  data.table::rbindlist(out)
}
