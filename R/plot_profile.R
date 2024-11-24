# plot_base_segment <- function(profile_data) {
#   ggplot(profile_data, aes(x = index, y = melatonin)) +
#     # Plot all points and lines for the full profile
#     geom_line(color = "blue", size = 1) +
#     geom_point(color = "blue", size = 2) +
#     # Highlight points belonging to the base segment
#     geom_point(
#       data = profile_data %>% filter(base == 1),
#       aes(x = index, y = melatonin),
#       color = "red", size = 3
#     ) +
#     labs(
#       title = "Melatonin Profile with Base Segment Highlighted",
#       x = "Index",
#       y = "Melatonin Concentration"
#     ) +
#     theme_minimal()
# }
# plot_base_segment <- function(profile_data) {
#   ggplot(profile_data, aes(x = index, y = melatonin)) +
#     # Plot all points and lines for the full profile
#     geom_line(color = "blue", size = 1) +
#     geom_point(color = "blue", size = 2) +
#     # Highlight points belonging to the base segment
#     geom_point(
#       data = profile_data %>% filter(base == 1),
#       aes(x = index, y = melatonin),
#       color = "red", size = 3
#     ) +
#     geom_line(
#       data = profile_data %>% filter(base == 1),
#       aes(x = index, y = melatonin),
#       color = "red", size = 1
#     ) +
#     labs(
#       title = "Melatonin Profile with Base Segment Highlighted",
#       x = "Index",
#       y = "Melatonin Concentration"
#     ) +
#     theme_minimal()
# }
#plot_base_segment <- function(profile_data) {
#  profile_data <- adjust_time_axis(profile_data)
#  ggplot(profile_data, aes(x = adjusted_time, y = melatonin)) +
#    # Plot all points and lines for the full profile
#    geom_line(color = "blue", size = 1) +
#    geom_point(color = "blue", size = 2) +
#    # Highlight points belonging to the base segment
#    geom_point(
#      data = profile_data %>% filter(base == 1),
#      aes(x = time, y = melatonin),
#      color = "red", size = 3
#    ) +
#    geom_line(
#      data = profile_data %>% filter(base == 1),
#      aes(x = time, y = melatonin),
#      color = "red", size = 1
#    ) +
#    labs(
#      title = "Melatonin Profile with Base Segment Highlighted",
#      x = "Time",
#      y = "Melatonin Concentration"
#    ) +
#    theme_minimal()
#}

# last best version
plot_base_segment <- function(profile_data) {
  ggplot2::ggplot(profile_data, ggplot2::aes(x = .data$datetime, y = .data$melatonin)) +
    # Plot a single continuous line for the full profile (mapped to "Full Profile")
    ggplot2::geom_line(
      ggplot2::aes(color = "Full Profile", group = 1),
      size = 1.25
    ) +
    # Overlay points for the base segment (mapped to "Base Segment")
    ggplot2::geom_line(
      data = dplyr::filter(profile_data, .data$base == 1),
      ggplot2::aes(x = .data$datetime, y = .data$melatonin, color = "Base Segment"),
      size = 1.25
    ) +
    # Format the x-axis to show only time
    ggplot2::scale_x_datetime(
      labels = scales::date_format("%H:%M"),
      date_breaks = "2 hours"
    ) +
    # Add plot labels
    ggplot2::labs(
      title = "Melatonin Profile with Base Segment Highlighted",
      x = "Time",
      y = "Melatonin Concentration"
    ) +
    ggplot2::theme_minimal() +
    # Show the legend on the right
    ggplot2::theme(legend.position = "right")
}

# # last best plot_profile version
# plot_profile <- function(profile_data, show_segments = c("full", "base", "ascending")) {
#   # Start the plot with the full profile line
#   plot <- ggplot2::ggplot(profile_data, ggplot2::aes(x = .data$datetime, y = .data$melatonin)) +
#     ggplot2::geom_line(
#       ggplot2::aes(color = "Full Profile", group = 1),
#       size = 1.25
#     ) +
#     ggplot2::scale_x_datetime(
#       labels = scales::date_format("%H:%M"),
#       date_breaks = "2 hours"
#     ) +
#     ggplot2::labs(
#       title = "Melatonin Profile with Segments",
#       x = "Time",
#       y = "Melatonin Concentration"
#     ) +
#     ggplot2::theme_minimal() +
#     ggplot2::theme(legend.position = "right")
#
#   # Add base segment line if required
#   if ("base" %in% show_segments) {
#     plot <- plot + ggplot2::geom_line(
#       data = dplyr::filter(profile_data, .data$base == 1),
#       ggplot2::aes(x = .data$datetime, y = .data$melatonin, color = "Base Segment"),
#       size = 1.25
#     )
#   }
#
#   # Add ascending segment line if required
#   if ("ascending" %in% show_segments) {
#     plot <- plot + ggplot2::geom_line(
#       data = dplyr::filter(profile_data, .data$ascending == 1),
#       ggplot2::aes(x = .data$datetime, y = .data$melatonin, color = "Ascending Segment"),
#       size = 1.25
#     )
#   }
#
#   # Return the plot
#   return(plot)
# }
#

plot_profile <- function(profile_data, show_segments = TRUE) {
  plot <- ggplot2::ggplot(profile_data, ggplot2::aes(x = .data$datetime, y = .data$melatonin)) +
    # Plot a single dotted line for the full profile (mapped to "Full Profile")
    ggplot2::geom_line(
      ggplot2::aes(color = "Full Profile", group = 1, linetype = "Full Profile"),
      size = 1.25
    ) +
    # Format the x-axis to show only time
    ggplot2::scale_x_datetime(
      labels = scales::date_format("%H:%M"),
      date_breaks = "2 hours"
    ) +
    # Add plot labels
    ggplot2::labs(
      title = "Melatonin Profile",
      x = "Time",
      y = "Melatonin Concentration"
    ) +
    ggplot2::theme_minimal() +
    # Show the legend on the right
    ggplot2::theme(legend.position = "right") +
    # Customize line types in the legend
    ggplot2::scale_linetype_manual(values = c("Full Profile" = "dotted"))

  # Add base and ascending segments if show_segments is TRUE
  if (show_segments) {
    plot <- plot +
      # Overlay points for the base segment (mapped to "Base Segment")
      ggplot2::geom_line(
        data = dplyr::filter(profile_data, .data$base == 1),
        ggplot2::aes(x = .data$datetime, y = .data$melatonin, color = "Base Segment"),
        size = 1.25
      ) +
      # Overlay points for the ascending segment (mapped to "Ascending Segment")
      ggplot2::geom_line(
        data = dplyr::filter(profile_data, .data$ascending == 1),
        ggplot2::aes(x = .data$datetime, y = .data$melatonin, color = "Ascending Segment"),
        size = 1.25
      )
  }

  # Return the plot object
  return(plot)
}


