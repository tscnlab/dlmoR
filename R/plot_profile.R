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
plot_base_segment <- function(profile_data) {
  ggplot2::ggplot(profile_data, ggplot2::aes(x = .data$datetime, y = .data$melatonin)) +
    # Plot a single continuous line for the full profile
    ggplot2::geom_line(ggplot2::aes(group = 1), color = "blue", size = 1.25) +
    # Overlay points for the base segment
    ggplot2::geom_line(
      data = dplyr::filter(profile_data, .data$base == 1),
      ggplot2::aes(x = .data$datetime, y = .data$melatonin),
      color = "red", size = 1.25
    ) +
    # Add plot labels and format the x-axis
    ggplot2::scale_x_datetime(
      labels = scales::date_format("%H:%M"),
      date_breaks = "2 hours"
    ) +
    ggplot2::labs(
      title = "Melatonin Profile with Base Segment Highlighted",
      x = "Time",
      y = "Melatonin Concentration"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")  # Optional: remove the legend
}




