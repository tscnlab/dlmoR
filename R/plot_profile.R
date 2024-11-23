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
plot_base_segment <- function(profile_data) {
  ggplot(profile_data, aes(x = datetime, y = melatonin)) +
    # Plot all points and lines for the full profile
    geom_line(color = "blue", size = 1) +
    geom_point(color = "blue", size = 2) +
    # Highlight points belonging to the base segment
    geom_point(
      data = profile_data %>% filter(base == 1),
      aes(x = datetime, y = melatonin),
      color = "red", size = 3
    ) +
    geom_line(
      data = profile_data %>% filter(base == 1),
      aes(x = datetime, y = melatonin),
      color = "red", size = 1
    ) +
    labs(
      title = "Melatonin Profile with Base Segment Highlighted",
      x = "Datetime",
      y = "Melatonin Concentration"
    ) +
    theme_minimal()
}
