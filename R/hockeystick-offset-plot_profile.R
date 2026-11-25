#' Plot a Hockey-Stick DLMOoff Result
#'
#' Draws a DLMOoff result in the same format as the hockey-stick DLMO plots
#' (the figure is made by [plot_profile()]): the whole input profile, the
#' threshold, base points (blue) and descending points (green; the "ascending"
#' segment of the mirrored fit), the fitted base line and descending line or
#' parabola, the DLMOoff (diamond), the region of interest with the residual
#' heatmap of the search, and the parallelogram.
#'
#' The title gives the DLMOoff clock time and date; the subtitle the fitted
#' functions `f_desc` and `f_base`, with `x_t` in decimal hours since 00:00 of
#' the DLMOoff date. On inputs longer than 26 h the time axis has 4-h ticks
#' with the date instead of 2-h ticks.
#'
#' @param result An object returned by
#'   `calculate_dlmo(method = "hockeystick", dlmo_off = TRUE)`.
#' @param plot_coarse Logical. Plot the coarse search (`TRUE`) or the fine
#'   search (`FALSE`, default; falls back to the coarse search if there is no
#'   fine search).
#'
#' @return A `ggplot2` object.
#'
#' @examples
#' filename <- system.file("extdata/sample_melatonin_profile.csv", package = "dlmoR")
#' off <- calculate_dlmo(file_path = filename, dlmo_off = TRUE, threshold = 2.3,
#'                       fine_flag = FALSE)
#' plot_profile_dlmoff(off, plot_coarse = TRUE)
#'
#' @seealso [calculate_dlmo()] (`method = "hockeystick"`, `dlmo_off = TRUE`)
#' @export
plot_profile_dlmoff <- function(result, plot_coarse = FALSE) {
  if (!is.list(result) || is.null(result$dlmoff) || is.null(result$ip)) {
    stop("`result` must come from calculate_dlmo(method = \"hockeystick\", dlmo_off = TRUE).",
         call. = FALSE)
  }
  lv <- if (plot_coarse || is.null(result$dlmoff$fine)) "coarse" else "fine"
  s <- lv
  ip <- result$ip
  grid_name <- c(coarse = "grid_big", fine = "grid_small")[[lv]]
  res_name <- c(coarse = "res_big", fine = "res_small")[[lv]]
  d <- result$dlmoff[[lv]]

  # plot_profile() works with the onset names: descending -> "ascending"
  fit_obj <- list(inflection_point = ip[[paste0("inflection_point_", s)]],
                  base_params = ip[[paste0("base_params_", s)]],
                  ascending_params = ip[[paste0("descending_params_", s)]],
                  grid = ip[[grid_name]], res = ip[[res_name]],
                  datetime_ref = ip$datetime_ref)
  d_onset <- list(time = d$time, fit_melatonin = d$fit_melatonin,
                  fit_lines = list(base = d$fit_lines$base, ascending = d$fit_lines$descending))
  dlmo_obj <- if (lv == "fine") list(coarse = NULL, fine = d_onset) else list(coarse = d_onset, fine = NULL)

  prof <- result$prof
  prof$ascending <- prof$descending
  # newest sample first: plot_fit() draws the base line from the first base row
  # to the inflection point and the descending fit to the last "ascending" row
  prof <- prof[order(prof$datetime, decreasing = TRUE), , drop = FALSE]

  p <- plot_profile(prof, show_threshold = TRUE, threshold = result$threshold,
                    show_segments = TRUE, show_parallelogram = !is.null(result$prl),
                    pll_result = result$prl, show_roi = TRUE, roi = result$roi,
                    show_dlmoIP = TRUE, dlmoFit = fit_obj, dlmo = dlmo_obj, show_fit = TRUE,
                    show_roi_heatmap = TRUE, plot_coarse = (lv == "coarse"))

  span_h <- as.numeric(difftime(max(prof$datetime), min(prof$datetime), units = "hours"))
  if (span_h > 26) {
    p <- suppressMessages(p + ggplot2::scale_x_datetime(
      labels = scales::date_format("%H:%M\n%d.%m"), date_breaks = "4 hours"))
  }

  fl <- d$fit_lines
  num <- function(x) as.numeric(x)
  subtitle <- if (identical(fl$descending$type, "parabolic")) {
    bquote(bold(f[desc]) == .(num(fl$descending$a)) * x[t]^2 + .(num(fl$descending$b)) * x[t] +
             .(num(fl$descending$c)) * "," ~ ~ ~ ~
             bold(f[base]) == .(num(fl$base$m)) * x[t] + .(num(fl$base$b)))
  } else {
    bquote(bold(f[desc]) == .(num(fl$descending$m)) * x[t] + .(num(fl$descending$b)) * "," ~ ~ ~ ~
             bold(f[base]) == .(num(fl$base$m)) * x[t] + .(num(fl$base$b)))
  }
  p + ggplot2::labs(
    title = paste0(if (lv == "fine") "Fine" else "Coarse", " Fit - DLMOff time: ",
                   format(d$datetime, "%H:%M:%S"), "  (", format(d$datetime, "%d.%m.%Y"), ")"),
    subtitle = subtitle)
}
