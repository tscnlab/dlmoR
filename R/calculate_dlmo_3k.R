#' Calculate DLMO with the 3k Threshold Method
#'
#' Estimates the dim-light melatonin onset (DLMO) with the 3k method: an
#' individual threshold is set at the mean of the first `n_baseline` samples
#' plus `k_sd` standard deviations (Voultsios et al., 1997), and the DLMO is
#' the time at which melatonin rises above this threshold and stays above it
#' for at least `min_duration` hours, located by linear interpolation
#' (Molina & Burgess, 2011).
#'
#' @param data A data frame with columns `datetime` (POSIXct) and `melatonin`
#'   (numeric), or `time` (hms) and `melatonin`. If `NULL`, `file_path` is used.
#' @param file_path Path to a CSV file with the same columns (comma or semicolon
#'   separated).
#' @param n_baseline Integer (>= 2). Number of initial samples that form the
#'   baseline (default 3).
#' @param k_sd Numeric (>= 0). Number of baseline standard deviations added to
#'   the baseline mean (default 2).
#' @param min_duration Numeric (>= 0). Hours melatonin must stay above the
#'   threshold for the crossing to count as the DLMO (default 2). Use 0 for the
#'   first crossing.
#' @param min_threshold Optional numeric. Lower bound for the threshold: if
#'   mean + `k_sd` x SD is below `min_threshold`, `min_threshold` is used
#'   instead (for example the assay limit of quantification). Useful when the
#'   baseline sits at the assay floor and its SD is close to 0. Default `NULL`
#'   (no lower bound, the original 3k method).
#' @param lod Optional numeric. Assay limit of detection or quantification;
#'   baseline values at or below it are reported as being at the assay floor.
#'   Used for warnings only.
#' @param max_gap Numeric (hours). A sampling gap longer than this within one
#'   hour before to `min_duration` hours after the DLMO triggers a warning
#'   (default 2).
#' @param plot Logical. If `TRUE` (default), a diagnostic plot is added to the
#'   result as `dlmoplot`.
#'
#' @return An object of class `dlmo_3k`, a list with
#' - **`dlmo_time`** (hms) and **`dlmo_datetime`** (POSIXct): the DLMO;
#' - **`status`**: `"ok"` (see the partial result attached to errors below);
#' - **`threshold`**: threshold used (pg/mL); **`threshold_raw`**: mean + k_sd x SD
#'   before `min_threshold`;
#' - **`baseline`**: list with the baseline samples (`datetime`, `melatonin`),
#'   their `mean`, `sd` and `n`;
#' - **`intervals`**: tibble of all intervals above the threshold (`start`, `end`,
#'   `duration_h`, `open_end`, ...); **`selected_interval`**: the row giving the DLMO;
#' - **`prof`**: the analysed profile with `hours_from_start`, `baseline` (0/1)
#'   and `above_threshold` (0/1) columns;
#' - **`issues`**: data frame of the problems found (`code`, `severity`,
#'   `message`); **`flags`**: named logical vector, one entry per issue code;
#' - **`settings`**: the arguments used;
#' - **`dlmoplot`**: `ggplot` object (if `plot = TRUE`).
#'
#' @details
#' **Errors.** The function stops with an error when
#' - the input is invalid: no input, invalid arguments, negative
#'   concentrations, duplicated time stamps or fewer than `n_baseline + 2`
#'   samples (see [check_3k_input()]);
#' - no DLMO can be determined (see [select_dlmo_3k()]):
#'   - `no_crossing`: melatonin never rises above the threshold;
#'   - `not_sustained`: it rises above the threshold but never stays above it
#'     for `min_duration` hours;
#'   - `truncated`: it rises above the threshold but the data end before
#'     `min_duration` hours.
#'
#' The "no DLMO" errors have class `dlmo_3k_error` and carry the partial result
#' as `result` (same elements as above, with `NA` DLMO, the failure `status`
#' and, if `plot = TRUE`, the plot), so the profile can still be inspected; see
#' the examples.
#'
#' **Warnings** (DLMO returned, but check it):
#' `baseline_sd_zero`, `baseline_at_floor` (see [check_3k_baseline()]),
#' `gap_near_dlmo` (long sampling gap around the DLMO) and
#' `missing_values_removed`.
#'
#' **Notes** (messages): `unsorted_input`, `threshold_floored` (`min_threshold`
#' was used) and `transient_excursion` (shorter excursions above the threshold
#' before the DLMO were ignored).
#'
#' A sample counts as above the threshold only when it is strictly greater than
#' the threshold, so a baseline with SD 0 (threshold equal to the baseline
#' value) is not mistaken for the rise; see [find_threshold_intervals()].
#'
#' Because the threshold comes from the first samples, the profile should start
#' in the daytime trough, a few hours before the expected onset. Profiles that
#' start during the morning decline or already on the rise give a biased
#' threshold; cut such profiles before calling this function.
#'
#' @references
#' Voultsios A, Kennaway DJ, Dawson D (1997). Salivary melatonin as a circadian
#' phase marker: validation and comparison to plasma melatonin. *J Biol Rhythms*
#' 12(5):457-466.
#'
#' Molina TA, Burgess HJ (2011). Calculating the dim light melatonin onset: the
#' impact of threshold and sampling rate. *Chronobiol Int* 28(8):714-718.
#'
#' @examples
#' filename <- system.file("extdata/sample_melatonin_profile.csv", package = "dlmoR")
#' sample_profile <- readr::read_delim(filename, delim = ";", show_col_types = FALSE)
#'
#' # The sample profile starts at 07:30, during the morning decline. Start it in
#' # the afternoon trough so that the baseline is taken there:
#' afternoon <- sample_profile[sample_profile$datetime >=
#'                               as.POSIXct("2024-04-17 12:00:00", tz = "UTC"), ]
#' res <- calculate_dlmo_3k(data = afternoon)
#' res
#' res$dlmoplot
#'
#' # With the morning samples the threshold is inflated and no DLMO is found.
#' # The error carries the partial result, so the profile can still be plotted:
#' failed <- tryCatch(calculate_dlmo_3k(data = sample_profile),
#'                    dlmo_3k_error = function(e) e$result)
#' failed$status
#' failed$dlmoplot
#'
#' @seealso [plot_profile_3k()], [compute_3k_threshold()],
#'   [find_threshold_intervals()], [select_dlmo_3k()], [check_3k_baseline()],
#'   [calculate_dlmo()] for the hockey-stick method.
#' @export
calculate_dlmo_3k <- function(data = NULL, file_path = NULL, n_baseline = 3, k_sd = 2,
                              min_duration = 2, min_threshold = NULL, lod = NULL,
                              max_gap = 2, plot = TRUE) {

  # ---- input -------------------------------------------------------------------
  if (is.null(data) && is.null(file_path)) {
    stop("You must provide either `data` or `file_path`.", call. = FALSE)
  }
  if (!is.null(file_path)) {
    message("Loading data from file: ", file_path)
    data <- .read_melatonin_data(file_path)
  }
  data <- validate_df_structure(data)
  data <- .ensure_datetime_3k(data)

  chk <- check_3k_input(data, n_baseline = n_baseline, k_sd = k_sd,
                        min_duration = min_duration, min_threshold = min_threshold,
                        lod = lod, max_gap = max_gap)
  prof <- chk$profile
  issues <- chk$issues
  n <- nrow(prof)
  t0 <- prof$datetime[1]
  t_h <- as.numeric(difftime(prof$datetime, t0, units = "hours"))
  to_datetime <- function(h) t0 + h * 3600

  # ---- threshold and baseline checks -----------------------------------------------
  thr <- compute_3k_threshold(prof$melatonin, n_baseline = n_baseline, k_sd = k_sd,
                              min_threshold = min_threshold)
  bidx <- thr$baseline_index
  base_issues <- check_3k_baseline(prof$melatonin[bidx], thr$threshold, lod = lod)
  if (thr$floored) {
    # the floor already addresses a threshold inside assay noise
    base_issues <- base_issues[!base_issues$code %in% c("baseline_sd_zero", "baseline_at_floor"), ]
    issues <- .add_issue(issues, "threshold_floored", sprintf(
      "Threshold raised from %.3f (mean + %g x SD) to `min_threshold` = %g.",
      thr$threshold_raw, k_sd, thr$threshold), severity = "note")
  }
  issues <- rbind(issues, base_issues)

  # ---- crossings and DLMO ------------------------------------------------------
  iv <- find_threshold_intervals(t_h, prof$melatonin, thr$threshold)
  sel <- select_dlmo_3k(iv, min_duration = min_duration)
  hhmm <- function(h) format(to_datetime(h), "%H:%M")

  failure <- switch(sel$status,
    no_crossing = sprintf(
      "No DLMO: melatonin never rises above the threshold (%.2f; profile maximum %.2f).",
      thr$threshold, max(prof$melatonin)),
    not_sustained = sprintf(
      paste0("No DLMO: melatonin rises above the threshold (%.2f) %d time(s) but never ",
             "stays above it for %g h (longest %.1f h)."),
      thr$threshold, nrow(iv), min_duration, max(iv$duration_h)),
    truncated = sprintf(
      paste0("No DLMO: melatonin rises above the threshold (%.2f) at %s, but the profile ",
             "ends %d min later, before the required %g h. Extend the profile or lower ",
             "`min_duration`."),
      thr$threshold, hhmm(iv$start_h[nrow(iv)]),
      as.integer(floor(iv$duration_h[nrow(iv)] * 60)), min_duration),
    NULL)

  if (sel$status == "ok") {
    if (sel$n_transient > 0) {
      tr <- iv[seq_len(sel$n_transient), ]
      issues <- .add_issue(issues, "transient_excursion", sprintf(
        "%d shorter excursion(s) above the threshold before the DLMO were ignored (starting %s).",
        sel$n_transient, paste(hhmm(tr$start_h), collapse = ", ")), severity = "note")
    }
    d <- sel$dlmo_h
    lo <- max(c(1, which(t_h <= d - 1)))
    hi <- which(t_h >= d + min_duration)
    hi <- if (length(hi)) hi[1] else n
    gaps <- diff(t_h[lo:hi])
    if (length(gaps) && max(gaps) > max_gap) {
      g <- which.max(gaps) + lo - 1
      issues <- .add_issue(issues, "gap_near_dlmo", sprintf(
        "Sampling gap of %.1f h near the DLMO (%s to %s); the interpolated DLMO is uncertain.",
        max(gaps), hhmm(t_h[g]), hhmm(t_h[g + 1])))
    }
  }

  # ---- report issues -----------------------------------------------------------
  for (i in seq_len(nrow(issues))) {
    if (issues$severity[i] == "warning") {
      warning("[3k] ", issues$message[i], call. = FALSE)
    } else {
      message("[3k] ", issues$message[i])
    }
  }
  codes <- c("missing_values_removed", "unsorted_input", "threshold_floored",
             "baseline_sd_zero", "baseline_at_floor", "transient_excursion",
             "gap_near_dlmo")
  flags <- stats::setNames(codes %in% issues$code, codes)

  # ---- result ------------------------------------------------------------------
  prof$hours_from_start <- t_h
  prof$baseline <- as.integer(seq_len(n) %in% bidx)
  prof$above_threshold <- as.integer(prof$melatonin > thr$threshold)

  intervals <- tidyr::as_tibble(data.frame(
    start = to_datetime(iv$start_h), end = to_datetime(iv$end_h),
    duration_h = iv$duration_h, start_interpolated = iv$start_interpolated,
    end_interpolated = iv$end_interpolated, open_end = iv$open_end,
    selected = seq_len(nrow(iv)) %in% sel$interval))

  dlmo_datetime <- to_datetime(sel$dlmo_h)

  result <- structure(list(
    dlmo_time = hms::as_hms(dlmo_datetime),
    dlmo_datetime = dlmo_datetime,
    status = sel$status,
    threshold = thr$threshold,
    threshold_raw = thr$threshold_raw,
    baseline = list(datetime = prof$datetime[bidx], melatonin = prof$melatonin[bidx],
                    mean = thr$baseline_mean, sd = thr$baseline_sd, n = n_baseline),
    intervals = intervals,
    selected_interval = sel$interval,
    prof = prof,
    issues = issues,
    flags = flags,
    settings = list(n_baseline = n_baseline, k_sd = k_sd, min_duration = min_duration,
                    min_threshold = min_threshold, lod = lod, max_gap = max_gap)
  ), class = c("dlmo_3k", "list"))

  if (plot) result$dlmoplot <- plot_profile_3k(result)

  if (!is.null(failure)) {
    stop(structure(
      class = c("dlmo_3k_error", "error", "condition"),
      list(message = paste0(failure, " The partial result is attached to the error ",
                            "(see ?calculate_dlmo_3k)."),
           call = NULL, result = result)))
  }
  result
}


#' Print a 3k DLMO Result
#'
#' @param x An object returned by [calculate_dlmo_3k()].
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.dlmo_3k <- function(x, ...) {
  cat("3k DLMO estimate\n")
  if (x$status == "ok") {
    cat("  DLMO:      ", format(x$dlmo_datetime, "%Y-%m-%d %H:%M:%S"), "\n")
  } else {
    cat("  DLMO:       not determined (", x$status, ")\n", sep = "")
  }
  cat(sprintf("  Threshold:  %.3f (baseline mean %.3f + %g x SD %.3f, n = %d)%s\n",
              x$threshold, x$baseline$mean, x$settings$k_sd, x$baseline$sd, x$baseline$n,
              if (isTRUE(x$flags[["threshold_floored"]])) " [raised to min_threshold]" else ""))
  cat(sprintf("  Duration:   above threshold for >= %g h\n", x$settings$min_duration))
  w <- x$issues$code[x$issues$severity == "warning"]
  if (length(w)) cat("  Warnings:  ", paste(w, collapse = ", "), "\n")
  invisible(x)
}
