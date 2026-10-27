# Load and prepare the base datasets into one common patient-level format:
# covariates, planned treatment time `s` (Inf if none), follow-up `time`, `status`,
# observed censoring time `cens_time`, administrative end `admin_time`, loss-to-follow-up
# indicator `cens_event` (outcome of the censoring model) and observed `trt_time`.

#' Prepare a base dataset
#' @param file Path to the raw file (read-only source).
#' @param spec Dataset entry of `cfg$datasets`.
#' @param name Data preparation key ("maringe" or "synpuf"; `spec$base` of the analysis).
#' @return data.table in the common format.
prepare_base_data <- function(file, spec, name) {
  checkmate::assert_file_exists(file)
  checkmate::assert_choice(name, c("maringe", "synpuf"))
  d <- switch(name, maringe = prepare_maringe(file), synpuf = prepare_synpuf(file, spec$horizon))
  # Treatment counts as received if it starts on or before the end of follow-up (tutorial rule)
  d[, trt_time := data.table::fifelse(s <= time, s, Inf)]
  d[]
}

#' Maringe et al. (2020) simulated lung cancer cohort
#' @param file CSV exported from ccwr::lungcancer.
#' @return data.table in the common format.
prepare_maringe <- function(file) {
  d <- data.table::fread(file)
  d[, .(
    id = as.character(id),
    age_c = (age - 80) / 5, # per 5 years around 80; keeps the Weibull fit well conditioned
    sex = as.integer(sex == 2),
    perf = factor(perf, levels = 0:2),
    stage = as.integer(stage == 2),
    deprivation,
    charlson = factor(charlson, levels = 0:2),
    emergency = as.integer(emergency),
    # Observed surgery date kept as the patient's (fixed) treatment plan; ADEMP D2
    s = data.table::fifelse(is.na(timetosurgery), Inf, as.numeric(timetosurgery)),
    time = as.numeric(fup_obs),
    status = as.integer(death),
    cens_time = 365.24,
    admin_time = 365.24,
    cens_event = 0L
  )]
}

#' Webster-Clark et al. (2025) SynPUF MI cohort (CCW tutorial base cohort)
#' @param file `baseccwcohort.sas7bdat` from the tutorial repository.
#' @param horizon Administrative end of follow-up (days).
#' @return data.table in the common format.
prepare_synpuf <- function(file, horizon) {
  d <- data.table::as.data.table(haven::read_sas(file))
  d[, .(
    id = ID,
    age_c = (baseline_age - 75) / 10,
    sex = as.integer(sex == "2"),
    renal = as.integer(renal == "Y"),
    nonwhite = as.integer(race != "1"),
    s = data.table::fifelse(is.na(RX_start), Inf, as.numeric(RX_start)),
    time = pmin(followup, horizon),
    status = as.integer(outcome == 1 & followup <= horizon),
    # Observed censoring: early loss to follow-up, otherwise the horizon
    cens_time = data.table::fifelse(outcome == 0 & followup < horizon, followup, horizon),
    admin_time = horizon,
    cens_event = as.integer(outcome == 0 & followup < horizon)
  )]
}
