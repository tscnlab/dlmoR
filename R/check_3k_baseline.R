#' Check the 3k Baseline for Known Problems
#'
#' The 3k threshold rests on only a few baseline samples, so its quality
#' depends on those samples sitting in a stable, low daytime trough. This
#' function flags baselines whose threshold is likely to lie within assay
#' noise. It never stops; it returns the problems found so that
#' [calculate_dlmo_3k()] can warn about them and record them.
#'
#' @param baseline_values Numeric vector of baseline melatonin values.
#' @param threshold The threshold used (after any `min_threshold`).
#' @param lod Optional assay limit of detection or quantification. If `NULL`,
#'   only exact zeros count as being at the assay floor.
#'
#' @return A data frame with columns `code`, `severity` (always `"warning"`)
#'   and `message`, one row per problem:
#'   - `baseline_sd_zero`: all baseline values are identical, so the threshold
#'     equals the baseline value and any increase counts as the rise;
#'   - `baseline_at_floor`: at least two baseline values are at the assay floor
#'     (0, or `<= lod`), so the threshold may lie within assay noise and the
#'     DLMO may be too early.
#'
#' @examples
#' check_3k_baseline(c(0, 0, 0), threshold = 0)
#' check_3k_baseline(c(0.2, 0.1, 0.3), threshold = 0.5, lod = 0.3)
#'
#' @seealso [calculate_dlmo_3k()]
#' @export
check_3k_baseline <- function(baseline_values, threshold, lod = NULL) {
  issues <- .new_issues()
  s <- stats::sd(baseline_values)
  vals <- paste(format(signif(baseline_values, 3), trim = TRUE), collapse = ", ")

  if (isTRUE(s == 0)) {
    issues <- .add_issue(issues, "baseline_sd_zero", sprintf(
      paste0("Baseline SD is 0 (baseline values: %s), so the threshold equals the baseline ",
             "value and any increase counts as the rise. Consider `min_threshold` ",
             "(e.g. the assay limit of quantification)."), vals))
  }

  floor_value <- if (is.null(lod)) 0 else lod
  n_floor <- sum(baseline_values <= floor_value)
  if (n_floor >= 2 && !isTRUE(s == 0)) {
    issues <- .add_issue(issues, "baseline_at_floor", sprintf(
      paste0("%d of %d baseline values are at the assay floor (<= %s; baseline values: %s). ",
             "The threshold (%.2f) may lie within assay noise, so the DLMO may be too early."),
      n_floor, length(baseline_values), format(floor_value), vals, threshold))
  }
  issues
}
