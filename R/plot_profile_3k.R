#' Plot a 3k DLMO Result
#'
#' Draws the melatonin profile with the baseline samples, the 3k threshold, the
#' intervals above the threshold and the DLMO.
#'
#' @param result An object returned by [calculate_dlmo_3k()].
#' @param show_intervals Logical. Shade the intervals above the threshold
#'   (default `TRUE`); the interval that gives the DLMO is shaded darker.
#'
#' @return A `ggplot2` object.
#'
#' @examples
#' filename <- system.file("extdata/sample_melatonin_profile.csv", package = "dlmoR")
#' sample_profile <- readr::read_delim(filename, delim = ";", show_col_types = FALSE)
#' afternoon <- sample_profile[sample_profile$datetime >=
#'                               as.POSIXct("2024-04-17 12:00:00", tz = "UTC"), ]
#' res <- calculate_dlmo_3k(data = afternoon, plot = FALSE)
#' plot_profile_3k(res)
#'
#' @seealso [calculate_dlmo_3k()]
#' @export
plot_profile_3k <- function(result, show_intervals = TRUE) {
  if (!inherits(result, "dlmo_3k")) stop("`result` must come from calculate_dlmo_3k().", call. = FALSE)
  prof <- result$prof
  base <- prof[prof$baseline == 1, ]

  title <- if (result$status == "ok") {
    paste("3k method - DLMO time:", format(result$dlmo_datetime, "%H:%M:%S"))
  } else {
    paste0("3k method - no DLMO (", result$status, ")")
  }
  subtitle <- sprintf("Threshold %.2f = baseline mean %.2f + %g x SD %.2f; above threshold for >= %g h",
                      result$threshold, result$baseline$mean, result$settings$k_sd,
                      result$baseline$sd, result$settings$min_duration)
  if (isTRUE(result$flags[["threshold_floored"]])) {
    subtitle <- sprintf("Threshold %.2f (min_threshold; raw %.2f); above threshold for >= %g h",
                        result$threshold, result$threshold_raw, result$settings$min_duration)
  }
  warn <- result$issues$code[result$issues$severity == "warning"]
  if (length(warn)) subtitle <- paste0(subtitle, "\nWarnings: ", paste(warn, collapse = ", "))

  p <- ggplot2::ggplot(prof, ggplot2::aes(x = .data$datetime, y = .data$melatonin))

  if (show_intervals && nrow(result$intervals)) {
    p <- p + ggplot2::geom_rect(
      data = result$intervals, inherit.aes = FALSE,
      ggplot2::aes(xmin = .data$start, xmax = .data$end, ymin = -Inf, ymax = Inf,
                   alpha = .data$selected),
      fill = "lightgreen") +
      ggplot2::scale_alpha_manual(values = c(`TRUE` = 0.35, `FALSE` = 0.15), guide = "none")
  }

  p <- p +
    ggplot2::geom_hline(yintercept = result$threshold, color = "burlywood3", linewidth = 1) +
    ggplot2::geom_line(color = "grey", linetype = "dotted", linewidth = 1) +
    ggplot2::geom_point(color = "grey40", size = 2) +
    ggplot2::geom_point(data = base, color = "#56B4E9", size = 3.5, shape = 21, stroke = 1.2)

  if (result$status == "ok") {
    p <- p +
      ggplot2::geom_vline(xintercept = result$dlmo_datetime, color = "deeppink", linewidth = 0.8) +
      ggplot2::annotate("point", x = result$dlmo_datetime, y = result$threshold,
                        color = "deeppink", size = 3)
  }

  tz <- attr(prof$datetime, "tzone")
  if (is.null(tz) || !nzchar(tz[1])) tz <- "UTC"
  p +
    ggplot2::scale_x_datetime(labels = scales::date_format("%H:%M", tz = tz[1]),
                              date_breaks = "2 hours") +
    ggplot2::labs(title = title, subtitle = subtitle, x = "Time", y = "Melatonin (pg/mL)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 9))
}
