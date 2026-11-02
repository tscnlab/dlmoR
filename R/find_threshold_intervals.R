#' Find Intervals Above a Threshold
#'
#' Finds every interval during which melatonin is above `threshold`. Start and
#' end times are found by linear interpolation between the samples on either
#' side of each crossing.
#'
#' @param time_h Numeric vector of sample times in hours (increasing, for
#'   example hours since the first sample).
#' @param melatonin Numeric vector of melatonin concentrations, same length.
#' @param threshold Single numeric threshold. A sample counts as above the
#'   threshold when it is strictly greater than `threshold`.
#'
#' @return A data frame with one row per interval:
#'   - `start_h`, `end_h`, `duration_h`: interval bounds and length in hours;
#'   - `start_interpolated`, `end_interpolated`: `FALSE` when the bound is a
#'     sample time rather than an interpolated crossing (profile starts or ends
#'     above the threshold);
#'   - `open_end`: `TRUE` if the profile is still above the threshold at its
#'     last sample, so the true end of the interval is unknown;
#'   - `first_index`, `last_index`: indices of the first and last sample above
#'     the threshold.
#'
#' @details
#' Samples equal to the threshold count as not above. This matters when the
#' baseline standard deviation is 0 (all baseline samples identical, often
#' values coded at the assay floor): the baseline samples then sit exactly on
#' the threshold and are not counted as a rise.
#'
#' @examples
#' find_threshold_intervals(0:5, c(1, 1, 3, 4, 1, 5), threshold = 2)
#'
#' @seealso [calculate_dlmo_3k()], [select_dlmo_3k()]
#' @export
find_threshold_intervals <- function(time_h, melatonin, threshold) {
  stopifnot(length(time_h) == length(melatonin), length(threshold) == 1)
  empty <- data.frame(start_h = numeric(0), end_h = numeric(0), duration_h = numeric(0),
                      start_interpolated = logical(0), end_interpolated = logical(0),
                      open_end = logical(0), first_index = integer(0), last_index = integer(0))
  above <- melatonin > threshold
  if (!any(above)) return(empty)

  n <- length(melatonin)
  cross <- function(i) {      # interpolated crossing between samples i and i + 1
    time_h[i] + (threshold - melatonin[i]) * (time_h[i + 1] - time_h[i]) /
      (melatonin[i + 1] - melatonin[i])
  }
  runs <- rle(above)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1

  out <- lapply(which(runs$values), function(j) {
    i0 <- starts[j]; i1 <- ends[j]
    s <- if (i0 > 1) cross(i0 - 1) else time_h[i0]
    e <- if (i1 < n) cross(i1) else time_h[i1]
    data.frame(start_h = s, end_h = e, duration_h = e - s,
               start_interpolated = i0 > 1, end_interpolated = i1 < n,
               open_end = i1 == n, first_index = i0, last_index = i1)
  })
  do.call(rbind, out)
}
