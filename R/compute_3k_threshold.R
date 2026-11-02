#' Compute the 3k Threshold
#'
#' Computes the individual threshold of the 3k method: the mean of the first
#' `n_baseline` samples plus `k_sd` times their standard deviation
#' (Voultsios et al., 1997; Molina & Burgess, 2011).
#'
#' @param melatonin Numeric vector of melatonin concentrations in time order.
#' @param n_baseline Integer (>= 2). Number of initial samples used as baseline.
#' @param k_sd Numeric (>= 0). Multiplier of the baseline standard deviation.
#' @param min_threshold Optional numeric. If the computed threshold is lower,
#'   `min_threshold` is used instead (for example the assay limit of
#'   quantification).
#'
#' @return A list with
#'   - `baseline_index`: indices of the baseline samples;
#'   - `baseline_mean`, `baseline_sd`: mean and sample standard deviation
#'     (denominator n - 1) of the baseline samples;
#'   - `threshold_raw`: `baseline_mean + k_sd * baseline_sd`;
#'   - `threshold`: the threshold actually used (raw value, or `min_threshold`
#'     if that is higher);
#'   - `floored`: `TRUE` if `min_threshold` replaced the raw threshold.
#'
#' @examples
#' compute_3k_threshold(c(1.2, 1.5, 1.3, 2.0, 6.5, 12))
#' compute_3k_threshold(c(0, 0, 0, 0.4, 3, 9), min_threshold = 0.9)
#'
#' @seealso [calculate_dlmo_3k()]
#' @export
compute_3k_threshold <- function(melatonin, n_baseline = 3, k_sd = 2, min_threshold = NULL) {
  if (length(melatonin) < n_baseline) {
    stop("The profile has fewer samples than `n_baseline`.", call. = FALSE)
  }
  idx <- seq_len(n_baseline)
  base <- melatonin[idx]
  base_mean <- mean(base)
  base_sd <- stats::sd(base)
  raw <- base_mean + k_sd * base_sd
  floored <- !is.null(min_threshold) && raw < min_threshold
  list(baseline_index = idx,
       baseline_mean = base_mean,
       baseline_sd = base_sd,
       threshold_raw = raw,
       threshold = if (floored) min_threshold else raw,
       floored = floored)
}
