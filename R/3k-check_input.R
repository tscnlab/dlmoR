#' Validate and Prepare a Melatonin Profile for the 3k Method
#'
#' Checks the arguments and the profile used by [calculate_dlmo_3k()] and
#' returns a clean, time-ordered profile. Problems that make an estimate
#' impossible stop with an error; problems that can be repaired (missing
#' values, unsorted rows) are repaired and reported as issues.
#'
#' @param data A data frame with a `datetime` (POSIXct) and a `melatonin`
#'   (numeric) column, as returned by [validate_df_structure()].
#' @param n_baseline Integer (>= 2). Number of initial samples used as baseline.
#' @param k_sd Numeric (>= 0). Number of baseline standard deviations added to
#'   the baseline mean.
#' @param min_duration Numeric (>= 0). Hours the profile must stay above the
#'   threshold.
#' @param min_threshold Optional numeric (>= 0). Lower bound for the threshold.
#' @param lod Optional numeric (>= 0). Assay limit of detection or
#'   quantification, used only for warnings.
#' @param max_gap Numeric (> 0). Sampling gap (hours) around the DLMO above
#'   which a warning is issued.
#'
#' @return A list with
#'   - `profile`: tibble with `datetime`, `time` and `melatonin`, sorted by time,
#'     without missing values;
#'   - `issues`: data frame with columns `code`, `severity` (`"warning"` or
#'     `"note"`) and `message` (possibly empty).
#'
#' @details
#' Errors (no estimate possible):
#' - invalid arguments (`n_baseline` < 2 or not a whole number, negative
#'   `k_sd`, `min_duration`, `min_threshold` or `lod`, non-positive `max_gap`);
#' - negative melatonin concentrations;
#' - duplicated time stamps;
#' - fewer than `n_baseline + 2` non-missing samples.
#'
#' Repaired and reported:
#' - rows with missing `datetime` or `melatonin` are removed
#'   (`missing_values_removed`);
#' - rows not in time order are sorted (`unsorted_input`).
#'
#' @seealso [calculate_dlmo_3k()]
#' @export
check_3k_input <- function(data, n_baseline = 3, k_sd = 2, min_duration = 2,
                           min_threshold = NULL, lod = NULL, max_gap = 2) {

  .check_number(n_baseline, "n_baseline", min = 2)
  if (n_baseline != round(n_baseline)) {
    stop("`n_baseline` must be a whole number (got ", n_baseline, ").", call. = FALSE)
  }
  .check_number(k_sd, "k_sd", min = 0)
  .check_number(min_duration, "min_duration", min = 0)
  if (!is.null(min_threshold)) .check_number(min_threshold, "min_threshold", min = 0)
  if (!is.null(lod)) .check_number(lod, "lod", min = 0)
  .check_number(max_gap, "max_gap", min = 0, strict = TRUE)

  if (!all(c("datetime", "melatonin") %in% names(data))) {
    stop("The profile must contain `datetime` and `melatonin` columns.", call. = FALSE)
  }

  issues <- .new_issues()

  keep <- !is.na(data$datetime) & !is.na(data$melatonin)
  if (any(!keep)) {
    issues <- .add_issue(issues, "missing_values_removed", sprintf(
      "%d row(s) with missing datetime or melatonin were removed.", sum(!keep)))
    data <- data[keep, , drop = FALSE]
  }

  if (any(data$melatonin < 0)) {
    stop("Melatonin concentrations must not be negative (found ",
         sum(data$melatonin < 0), " negative value(s)).", call. = FALSE)
  }

  if (anyDuplicated(data$datetime)) {
    dup <- unique(data$datetime[duplicated(data$datetime)])
    stop("Duplicated time stamps found: ",
         paste(format(dup, "%Y-%m-%d %H:%M:%S"), collapse = ", "),
         ". Remove or average duplicated samples first.", call. = FALSE)
  }

  if (is.unsorted(data$datetime)) {
    issues <- .add_issue(issues, "unsorted_input",
                         "Samples were not in time order and have been sorted.", severity = "note")
    data <- data[order(data$datetime), , drop = FALSE]
  }

  if (nrow(data) < n_baseline + 2) {
    stop(sprintf(paste0("At least %d non-missing samples are needed (%d baseline samples ",
                        "plus 2 to detect a crossing); found %d."),
                 n_baseline + 2, n_baseline, nrow(data)), call. = FALSE)
  }

  data$time <- hms::as_hms(data$datetime)
  list(profile = tidyr::as_tibble(data[, c("datetime", "time", "melatonin")]),
       issues = issues)
}


# ---- internal helpers shared by the 3k functions ------------------------------

.check_number <- function(x, name, min = -Inf, strict = FALSE) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x)) {
    stop("`", name, "` must be a single number.", call. = FALSE)
  }
  bad <- if (strict) x <= min else x < min
  if (bad) {
    stop("`", name, "` must be ", if (strict) "greater than " else "at least ",
         min, " (got ", x, ").", call. = FALSE)
  }
  invisible(TRUE)
}

.new_issues <- function() data.frame(code = character(0), severity = character(0),
                                     message = character(0), stringsAsFactors = FALSE)

.add_issue <- function(issues, code, message, severity = "warning") {
  rbind(issues, data.frame(code = code, severity = severity, message = message,
                           stringsAsFactors = FALSE))
}

# Ensure a POSIXct `datetime` column exists. When only clock times (hms) are
# given, a continuous time axis is built by adding a day at every midnight
# crossing (samples must be in collection order).
.ensure_datetime_3k <- function(data) {
  if ("datetime" %in% names(data)) return(data)
  secs <- as.numeric(data$time)
  day_offset <- cumsum(c(0, diff(secs) < 0)) * 86400
  message("No 'datetime' column: building one from 'time', assuming samples are in ",
          "collection order and adding a day at each midnight crossing.")
  data$datetime <- as.POSIXct("1970-01-01", tz = "UTC") + secs + day_offset
  data
}
