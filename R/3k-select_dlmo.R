#' Select the DLMO from Above-Threshold Intervals
#'
#' Applies the duration rule of the 3k method: the DLMO is the start of the
#' first interval that stays above the threshold for at least `min_duration`
#' hours (Molina & Burgess, 2011, used 2 hours).
#'
#' @param intervals Data frame returned by [find_threshold_intervals()].
#' @param min_duration Numeric (>= 0). Required time above the threshold, in hours.
#'
#' @return A list with
#'   - `dlmo_h`: start of the selected interval (same unit as the interval
#'     table), or `NA`;
#'   - `status`: one of
#'     - `"ok"`: a qualifying interval was found;
#'     - `"no_crossing"`: the profile never rises above the threshold;
#'     - `"not_sustained"`: it rises above the threshold but never for
#'       `min_duration` hours;
#'     - `"truncated"`: the last interval reaches the end of the data before
#'       `min_duration` hours, so it cannot be confirmed;
#'   - `interval`: row index of the selected interval, or `NA`;
#'   - `n_transient`: number of shorter excursions before the selected interval
#'     that were ignored.
#'
#' @examples
#' iv <- find_threshold_intervals(0:8, c(1, 1, 1, 4, 1, 3, 5, 6, 7), threshold = 2)
#' select_dlmo_3k(iv, min_duration = 2)
#'
#' @seealso [calculate_dlmo()] (`method = "3k"`), [find_threshold_intervals()]
#' @export
select_dlmo_3k <- function(intervals, min_duration = 2) {
  if (nrow(intervals) == 0) {
    return(list(dlmo_h = NA_real_, status = "no_crossing", interval = NA_integer_, n_transient = 0L))
  }
  ok <- which(intervals$duration_h >= min_duration)
  if (length(ok)) {
    k <- ok[1]
    return(list(dlmo_h = intervals$start_h[k], status = "ok", interval = k,
                n_transient = as.integer(k - 1)))
  }
  last <- nrow(intervals)
  status <- if (intervals$open_end[last]) "truncated" else "not_sustained"
  list(dlmo_h = NA_real_, status = status, interval = NA_integer_,
       n_transient = as.integer(nrow(intervals)))
}
