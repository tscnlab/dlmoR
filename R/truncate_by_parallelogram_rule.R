# # #
# # # # Helper function to plot the parallelogram and diagonals
# # # plot_parallelogram <- function(profile_data, params) {
# # #   # Ensure the datetime column is POSIXct (no change needed here since it's already in correct format)
# # #
# # #   # Get the ascending data
# # #   ascending_data <- profile_data %>% dplyr::filter(.data$ascending == 1)
# # #   x <- as.numeric(ascending_data$datetime)  # Convert datetime to numeric
# # #   y <- ascending_data$melatonin
# # #
# # #   # Calculate the parallelogram corners
# # #   y0 <- min(y)
# # #   y1 <- max(y)
# # #   corners <- get_corners(params[1], y0, params[2], y1, params[3])
# # #   lower_left <- corners[[1]]
# # #   lower_right <- corners[[2]]
# # #   upper_right <- corners[[3]]
# # #   upper_left <- corners[[4]]
# # #
# # #   # Create a tibble of the parallelogram corners
# # #   parallelogram_data <- tibble(
# # #     x = c(lower_left[1], lower_right[1], upper_right[1], upper_left[1], lower_left[1]),
# # #     y = c(lower_left[2], lower_right[2], upper_right[2], upper_left[2], lower_left[2])
# # #   )
# # #
# # #   # Create the base plot with the full profile line
# # #   plot <- ggplot2::ggplot(profile_data, ggplot2::aes(x = .data$datetime, y = .data$melatonin)) +
# # #     ggplot2::geom_line(
# # #       ggplot2::aes(color = "Full Profile", group = 1, linetype = "Full Profile"),
# # #       size = 1.25
# # #     ) +
# # #     ggplot2::scale_x_datetime(
# # #       labels = scales::date_format("%H:%M"),
# # #       date_breaks = "2 hours"
# # #     ) +
# # #     ggplot2::labs(
# # #       title = "Melatonin Profile with Parallelogram",
# # #       x = "Time",
# # #       y = "Melatonin Concentration"
# # #     ) +
# # #     ggplot2::theme_minimal() +
# # #     ggplot2::theme(legend.position = "right") +
# # #     ggplot2::scale_linetype_manual(values = c("Full Profile" = "dotted"))
# # #
# # #   # Add base and ascending segments
# # #   plot <- plot +
# # #     ggplot2::geom_line(
# # #       data = dplyr::filter(profile_data, .data$base == 1),
# # #       ggplot2::aes(x = .data$datetime, y = .data$melatonin, color = "Base Segment"),
# # #       size = 1.25
# # #     ) +
# # #     ggplot2::geom_line(
# # #       data = dplyr::filter(profile_data, .data$ascending == 1),
# # #       ggplot2::aes(x = .data$datetime, y = .data$melatonin, color = "Ascending Segment"),
# # #       size = 1.25
# # #     )
# # #
# # #   # Plot the parallelogram
# # #   plot <- plot +
# # #     ggplot2::geom_polygon(
# # #       data = parallelogram_data,
# # #       ggplot2::aes(x = .data$x, y = .data$y),
# # #       fill = "blue", alpha = 0.2, color = "black"
# # #     )
# # #
# # #   # Use annotate to plot the diagonals of the parallelogram
# # #   plot <- plot +
# # #     ggplot2::annotate("segment",
# # #                       x = lower_left[1], y = lower_left[2],
# # #                       xend = upper_right[1], yend = upper_right[2],
# # #                       color = "red", linetype = "dashed", size = 1
# # #     ) +
# # #     ggplot2::annotate("segment",
# # #                       x = lower_right[1], y = lower_right[2],
# # #                       xend = upper_left[1], yend = upper_left[2],
# # #                       color = "red", linetype = "dashed", size = 1
# # #     )
# # #
# # #   return(plot)
# # # }
# # #
# # # # Helper function to compute parallelogram corners
# # # get_corners <- function(x0, y0, x1, y1, slope) {
# # #   if (slope == 0) {
# # #     lower_left <- c(x0, y0)
# # #     lower_right <- c(x1, y0)
# # #     upper_left <- c(x0, y1)
# # #     upper_right <- c(x1, y1)
# # #   } else {
# # #     height <- y1 - y0
# # #     delta_x <- height / slope
# # #
# # #     if (abs(slope) > 1e3) {
# # #       delta_x <- 0
# # #     }
# # #
# # #     lower_left <- c(x0, y0)
# # #     lower_right <- c(x1, y0)
# # #     upper_left <- c(x0 + delta_x, y1)
# # #     upper_right <- c(x1 + delta_x, y1)
# # #   }
# # #
# # #   return(list(lower_left, lower_right, upper_right, upper_left))
# # # }
# # #
# # # # Helper function to compute constraints
# # # constraints <- function(params, x, y, y0, y1) {
# # #   x0 <- params[1]
# # #   x1 <- params[2]
# # #   slope <- params[3]
# # #   constraint_vals <- numeric()
# # #
# # #   corners <- get_corners(x0, y0, x1, y1, slope)
# # #   lower_left <- corners[[1]]
# # #   lower_right <- corners[[2]]
# # #   upper_right <- corners[[3]]
# # #   upper_left <- corners[[4]]
# # #
# # #   for (i in seq_along(x)) {
# # #     xi <- x[i]
# # #     yi <- y[i]
# # #
# # #     y_lower <- if (xi <= lower_right[1]) y0 else y0 + slope * (xi - lower_right[1])
# # #     y_upper <- if (xi >= upper_left[1]) y1 else y1 + slope * (xi - upper_left[1])
# # #
# # #     x_lower <- min(lower_left[1], upper_left[1])
# # #     x_upper <- max(upper_right[1], lower_right[1])
# # #
# # #     constraint_vals <- c(constraint_vals, y_upper - yi)
# # #     constraint_vals <- c(constraint_vals, yi - y_lower)
# # #     constraint_vals <- c(constraint_vals, xi - x_lower)
# # #     constraint_vals <- c(constraint_vals, x_upper - xi)
# # #   }
# # #   return(constraint_vals)
# # # }
# # #
# # # # Helper function for objective minimization
# # # objective <- function(params, x, y, y0, y1) {
# # #   x0 <- params[1]
# # #   x1 <- params[2]
# # #   slope <- params[3]
# # #   corners <- get_corners(x0, y0, x1, y1, slope)
# # #   v1 <- corners[[2]] - corners[[1]]
# # #   v2 <- corners[[4]] - corners[[1]]
# # #   area <- abs(v1[1] * v2[2] - v1[2] * v2[1])
# # #
# # #   c_penalty <- min(constraints(params, x, y, y0, y1))
# # #   if (c_penalty > 0) {
# # #     c_penalty <- 0
# # #   }
# # #   return(area + 1e3 * c_penalty^2)
# # # }
# # #
# # # # Function to optimize parallelogram
# # # optimize_parallelogram <- function(x, y) {
# # #   y0 <- min(y)
# # #   y1 <- max(y)
# # #   x0_initial <- min(x)
# # #   x1_initial <- max(x)
# # #   slope_initial <- (y1 - y0) / (x1_initial - x0_initial)
# # #
# # #   initial_guess <- c(x0_initial, x1_initial, slope_initial)
# # #
# # #   result <- optim(
# # #     par = initial_guess,
# # #     fn = objective,
# # #     method = "L-BFGS-B",
# # #     lower = c(-Inf, -Inf, 0),
# # #     upper = c(Inf, Inf, Inf),
# # #     x = x, y = y, y0 = y0, y1 = y1
# # #   )
# # #
# # #   return(result$par)
# # # }
# # #
# # # # Function to calculate diagonal slopes
# # # calculate_diagonal_slopes <- function(params, y0, y1) {
# # #   x0 <- params[1]
# # #   x1 <- params[2]
# # #   slope <- params[3]
# # #
# # #   corners <- get_corners(x0, y0, x1, y1, slope)
# # #   lower_left <- corners[[1]]
# # #   upper_right <- corners[[3]]
# # #   upper_left <- corners[[4]]
# # #   lower_right <- corners[[2]]
# # #
# # #   # Long diagonal
# # #   long_diag_slope <- (upper_right[2] - lower_left[2]) / (upper_right[1] - lower_left[1])
# # #
# # #   # Short diagonal
# # #   short_diag_slope <- (lower_right[2] - upper_left[2]) / (lower_right[1] - upper_left[1])
# # #
# # #   return(c(long_diag_slope = long_diag_slope, short_diag_slope = short_diag_slope))
# # # }
# # #
# # # # Main function to truncate by parallelogram rule
# # # truncate_by_parallelogram_rule <- function(profile_data) {
# # #   # Extract ascending segment
# # #   ascending_data <- profile_data %>% dplyr::filter(.data$ascending == 1)
# # #
# # #   if (nrow(ascending_data) < 2) {
# # #     warning("Insufficient ascending segment to apply parallelogram rule.")
# # #     return(profile_data)
# # #   }
# # #
# # #   x <- as.numeric(ascending_data$datetime)  # Convert datetime to numeric
# # #   y <- ascending_data$melatonin
# # #
# # #   # Optimize the parallelogram
# # #   params <- optimize_parallelogram(x, y)
# # # print(params)
# # #   # Calculate diagonal slopes
# # #   y0 <- min(y)
# # #   y1 <- max(y)
# # #   diagonals <- calculate_diagonal_slopes(params, y0, y1)
# # #   long_diag_slope <- diagonals["long_diag_slope"]
# # #   short_diag_slope <- diagonals["short_diag_slope"]
# # #
# # #   # Check the slope condition
# # #   if (abs(long_diag_slope) < 0.5 * abs(short_diag_slope)) {
# # #     # Truncate the ascending segment
# # #     profile_data <- profile_data %>%
# # #       dplyr::mutate(
# # #         ascending = ifelse(.data$ascending == 1 & .data$datetime > ascending_data$datetime[1], 0, .data$ascending)
# # #       )
# # #   }
# # #
# # #   # Plot the profile with parallelogram and diagonals
# # #   # plot <- plot_parallelogram(profile_data, params)
# # #   # print(plot)
# # #
# # #
# # #   return(profile_data)
# # # }
# # #
# # #
# # #
# # #
# # # # optimization works here
# # # # Helper function to plot the parallelogram and diagonals
# # # plot_parallelogram <- function(profile_data, params) {
# # #   # Ensure the datetime column is POSIXct (no change needed here since it's already in correct format)
# # #
# # #   # Get the ascending data
# # #   ascending_data <- profile_data %>% dplyr::filter(.data$ascending == 1)
# # #   x <- as.numeric(ascending_data$datetime)  # Convert datetime to numeric
# # #   y <- ascending_data$melatonin
# # #
# # #   # Calculate the parallelogram corners
# # #   y0 <- min(y)
# # #   y1 <- max(y)
# # #   corners <- get_corners(params[1], y0, params[2], y1, params[3])
# # #   lower_left <- corners[[1]]
# # #   lower_right <- corners[[2]]
# # #   upper_right <- corners[[3]]
# # #   upper_left <- corners[[4]]
# # #
# # #   # Create a tibble of the parallelogram corners
# # #   parallelogram_data <- tibble(
# # #     x = c(lower_left[1], lower_right[1], upper_right[1], upper_left[1], lower_left[1]),
# # #     y = c(lower_left[2], lower_right[2], upper_right[2], upper_left[2], lower_left[2])
# # #   )
# # #
# # #   # Create the base plot with the full profile line
# # #   plot <- ggplot2::ggplot(profile_data, ggplot2::aes(x = .data$datetime, y = .data$melatonin)) +
# # #     ggplot2::geom_line(
# # #       ggplot2::aes(color = "Full Profile", group = 1, linetype = "Full Profile"),
# # #       size = 1.25
# # #     ) +
# # #     ggplot2::scale_x_datetime(
# # #       labels = scales::date_format("%H:%M"),
# # #       date_breaks = "2 hours"
# # #     ) +
# # #     ggplot2::labs(
# # #       title = "Melatonin Profile with Parallelogram",
# # #       x = "Time",
# # #       y = "Melatonin Concentration"
# # #     ) +
# # #     ggplot2::theme_minimal() +
# # #     ggplot2::theme(legend.position = "right") +
# # #     ggplot2::scale_linetype_manual(values = c("Full Profile" = "dotted"))
# # #
# # #   # Add base and ascending segments
# # #   plot <- plot +
# # #     ggplot2::geom_line(
# # #       data = dplyr::filter(profile_data, .data$base == 1),
# # #       ggplot2::aes(x = .data$datetime, y = .data$melatonin, color = "Base Segment"),
# # #       size = 1.25
# # #     ) +
# # #     ggplot2::geom_line(
# # #       data = dplyr::filter(profile_data, .data$ascending == 1),
# # #       ggplot2::aes(x = .data$datetime, y = .data$melatonin, color = "Ascending Segment"),
# # #       size = 1.25
# # #     )
# # #
# # #   # Plot the parallelogram
# # #   plot <- plot +
# # #     ggplot2::geom_polygon(
# # #       data = parallelogram_data,
# # #       ggplot2::aes(x = .data$x, y = .data$y),
# # #       fill = "blue", alpha = 0.2, color = "black"
# # #     )
# # #
# # #   # Use annotate to plot the diagonals of the parallelogram
# # #   plot <- plot +
# # #     ggplot2::annotate("segment",
# # #                       x = lower_left[1], y = lower_left[2],
# # #                       xend = upper_right[1], yend = upper_right[2],
# # #                       color = "red", linetype = "dashed", size = 1
# # #     ) +
# # #     ggplot2::annotate("segment",
# # #                       x = lower_right[1], y = lower_right[2],
# # #                       xend = upper_left[1], yend = upper_left[2],
# # #                       color = "red", linetype = "dashed", size = 1
# # #     )
# # #
# # #   return(plot)
# # # }
# # #
# # # # Helper function to compute parallelogram corners
# # # get_corners <- function(x0, y0, x1, y1, slope) {
# # #   if (slope == 0) {
# # #     lower_left <- c(x0, y0)
# # #     lower_right <- c(x1, y0)
# # #     upper_left <- c(x0, y1)
# # #     upper_right <- c(x1, y1)
# # #   } else {
# # #     height <- y1 - y0
# # #     delta_x <- height / slope
# # #
# # #     if (abs(slope) > 1e3) {
# # #       delta_x <- 0
# # #     }
# # #
# # #     lower_left <- c(x0, y0)
# # #     lower_right <- c(x1, y0)
# # #     upper_left <- c(x0 + delta_x, y1)
# # #     upper_right <- c(x1 + delta_x, y1)
# # #   }
# # #
# # #   return(list(lower_left, lower_right, upper_right, upper_left))
# # # }
# # #
# # # # Helper function to compute constraints
# # # constraints <- function(params, x, y, y0, y1) {
# # #   x0 <- params[1]
# # #   x1 <- params[2]
# # #   slope <- params[3]
# # #   constraint_vals <- numeric()
# # #
# # #   corners <- get_corners(x0, y0, x1, y1, slope)
# # #   lower_left <- corners[[1]]
# # #   lower_right <- corners[[2]]
# # #   upper_right <- corners[[3]]
# # #   upper_left <- corners[[4]]
# # #
# # #   for (i in seq_along(x)) {
# # #     xi <- x[i]
# # #     yi <- y[i]
# # #
# # #     y_lower <- if (xi <= lower_right[1]) y0 else y0 + slope * (xi - lower_right[1])
# # #     y_upper <- if (xi >= upper_left[1]) y1 else y1 + slope * (xi - upper_left[1])
# # #
# # #     x_lower <- min(lower_left[1], upper_left[1])
# # #     x_upper <- max(upper_right[1], lower_right[1])
# # #
# # #     constraint_vals <- c(constraint_vals, y_upper - yi)
# # #     constraint_vals <- c(constraint_vals, yi - y_lower)
# # #     constraint_vals <- c(constraint_vals, xi - x_lower)
# # #     constraint_vals <- c(constraint_vals, x_upper - xi)
# # #   }
# # #   return(constraint_vals)
# # # }
# # #
# # # # Helper function for objective minimization
# # # objective <- function(params, x, y, y0, y1) {
# # #   x0 <- params[1]
# # #   x1 <- params[2]
# # #   slope <- params[3]
# # #   corners <- get_corners(x0, y0, x1, y1, slope)
# # #   v1 <- corners[[2]] - corners[[1]]
# # #   v2 <- corners[[4]] - corners[[1]]
# # #   area <- abs(v1[1] * v2[2] - v1[2] * v2[1])
# # #
# # #   c_penalty <- min(constraints(params, x, y, y0, y1))
# # #   if (c_penalty > 0) {
# # #     c_penalty <- 0
# # #   }
# # #   return(area + 1e3 * c_penalty^2)
# # # }
# # #
# # # # Function to optimize parallelogram using indices
# # # optimize_parallelogram <- function(profile_data) {
# # #   # Extract ascending segment and convert time to indices
# # #   ascending_data <- profile_data %>% dplyr::filter(.data$ascending == 1)
# # #
# # #   if (nrow(ascending_data) < 2) {
# # #     warning("Insufficient ascending segment to apply parallelogram rule.")
# # #     return(NULL)
# # #   }
# # #
# # #   x <- seq_along(ascending_data$melatonin)  # Use indices instead of datetime values
# # #   y <- ascending_data$melatonin
# # #
# # #   # Initial guess for x0, x1, and slope
# # #   y0 <- min(y)
# # #   y1 <- max(y)
# # #   x0_initial <- 1
# # #   x1_initial <- length(x)
# # #   slope_initial <- (y1 - y0) / (x1_initial - x0_initial)
# # #
# # #   initial_guess <- c(x0_initial, x1_initial, slope_initial)
# # #
# # #   result <- optim(
# # #     par = initial_guess,
# # #     fn = objective,
# # #     method = "L-BFGS-B",
# # #     lower = c(1, 1, -Inf),
# # #     upper = c(length(x), length(x), Inf),
# # #     x = x, y = y, y0 = y0, y1 = y1
# # #   )
# # #
# # #   return(result$par)
# # # }
# # #
# # # # Main function to truncate by parallelogram rule
# # # truncate_by_parallelogram_rule <- function(profile_data) {
# # #   # Extract ascending segment
# # #   ascending_data <- profile_data %>% dplyr::filter(.data$ascending == 1)
# # #
# # #   if (nrow(ascending_data) < 2) {
# # #     warning("Insufficient ascending segment to apply parallelogram rule.")
# # #     return(profile_data)
# # #   }
# # #
# # #   # Optimize the parallelogram
# # #   params <- optimize_parallelogram(profile_data)
# # #   if (is.null(params)) return(profile_data)
# # #   print(params)
# # #
# # #   # Calculate diagonal slopes
# # #   y0 <- min(ascending_data$melatonin)
# # #   y1 <- max(ascending_data$melatonin)
# # #   diagonals <- calculate_diagonal_slopes(params, y0, y1)
# # #   long_diag_slope <- diagonals["long_diag_slope"]
# # #   short_diag_slope <- diagonals["short_diag_slope"]
# # #
# # #   # Check the slope condition
# # #   if (abs(long_diag_slope) < 0.5 * abs(short_diag_slope)) {
# # #     # Truncate the ascending segment
# # #     profile_data <- profile_data %>%
# # #       dplyr::mutate(
# # #         ascending = ifelse(.data$ascending == 1 & .data$datetime > ascending_data$datetime[1], 0, .data$ascending)
# # #       )
# # #   }
# # #
# # #   # Plot the profile with parallelogram and diagonals
# # #   # plot <- plot_parallelogram(profile_data, params)
# # #   # print(plot)
# # #
# # #   return(profile_data)
# # # }
# # # optimization likely okay below:
# #
# # # Helper function to plot the parallelogram and diagonals
# # # plot_parallelogram <- function(profile_data, params) {
# # #   # Ensure the datetime column is POSIXct (no change needed here since it's already in correct format)
# # #
# # #   # Get the ascending data
# # #   ascending_data <- profile_data %>% dplyr::filter(.data$ascending == 1)
# # #   x <- as.numeric(ascending_data$datetime)  # Convert datetime to numeric
# # #   y <- ascending_data$melatonin
# # #
# # #   # Calculate the parallelogram corners
# # #   y0 <- min(y)
# # #   y1 <- max(y)
# # #   corners <- get_corners(params[1], y0, params[2], y1, params[3])
# # #   lower_left <- corners[[1]]
# # #   lower_right <- corners[[2]]
# # #   upper_right <- corners[[3]]
# # #   upper_left <- corners[[4]]
# # #
# # #   # Create a tibble of the parallelogram corners
# # #   parallelogram_data <- tibble(
# # #     x = c(lower_left[1], lower_right[1], upper_right[1], upper_left[1], lower_left[1]),
# # #     y = c(lower_left[2], lower_right[2], upper_right[2], upper_left[2], lower_left[2])
# # #   )
# # #
# # #   # Create the base plot with the full profile line
# # #   plot <- ggplot2::ggplot(profile_data, ggplot2::aes(x = .data$datetime, y = .data$melatonin)) +
# # #     ggplot2::geom_line(
# # #       ggplot2::aes(color = "Full Profile", group = 1, linetype = "Full Profile"),
# # #       size = 1.25
# # #     ) +
# # #     ggplot2::scale_x_datetime(
# # #       labels = scales::date_format("%H:%M"),
# # #       date_breaks = "2 hours"
# # #     ) +
# # #     ggplot2::labs(
# # #       title = "Melatonin Profile with Parallelogram",
# # #       x = "Time",
# # #       y = "Melatonin Concentration"
# # #     ) +
# # #     ggplot2::theme_minimal() +
# # #     ggplot2::theme(legend.position = "right") +
# # #     ggplot2::scale_linetype_manual(values = c("Full Profile" = "dotted"))
# # #
# # #   # Add base and ascending segments
# # #   plot <- plot +
# # #     ggplot2::geom_line(
# # #       data = dplyr::filter(profile_data, .data$base == 1),
# # #       ggplot2::aes(x = .data$datetime, y = .data$melatonin, color = "Base Segment"),
# # #       size = 1.25
# # #     ) +
# # #     ggplot2::geom_line(
# # #       data = dplyr::filter(profile_data, .data$ascending == 1),
# # #       ggplot2::aes(x = .data$datetime, y = .data$melatonin, color = "Ascending Segment"),
# # #       size = 1.25
# # #     )
# # #
# # #   # Plot the parallelogram
# # #   plot <- plot +
# # #     ggplot2::geom_polygon(
# # #       data = parallelogram_data,
# # #       ggplot2::aes(x = .data$x, y = .data$y),
# # #       fill = "blue", alpha = 0.2, color = "black"
# # #     )
# # #
# # #   # Use annotate to plot the diagonals of the parallelogram
# # #   plot <- plot +
# # #     ggplot2::annotate("segment",
# # #                       x = lower_left[1], y = lower_left[2],
# # #                       xend = upper_right[1], yend = upper_right[2],
# # #                       color = "red", linetype = "dashed", size = 1
# # #     ) +
# # #     ggplot2::annotate("segment",
# # #                       x = lower_right[1], y = lower_right[2],
# # #                       xend = upper_left[1], yend = upper_left[2],
# # #                       color = "red", linetype = "dashed", size = 1
# # #     )
# # #
# # #   return(plot)
# # # }
# # # Helper function to plot the parallelogram and diagonals
# # # plot_parallelogram <- function(profile_data, params) {
# # #   # Ensure the datetime column is POSIXct (no change needed here since it's already in correct format)
# # #
# # #   # Get the ascending data
# # #   ascending_data <- profile_data %>% dplyr::filter(.data$ascending == 1)
# # #   x <- as.numeric(ascending_data$datetime)  # Convert datetime to numeric
# # #   y <- ascending_data$melatonin
# # #
# # #   # Calculate the parallelogram corners (with datetime reference)
# # #   y0 <- min(y)
# # #   y1 <- max(y)
# # #   corners <- get_corners(params[1], y0, params[2], y1, params[3], ascending_data$datetime)
# # #   lower_left <- corners[[1]]
# # #   lower_right <- corners[[2]]
# # #   upper_right <- corners[[3]]
# # #   upper_left <- corners[[4]]
# # #
# # #   # Create a tibble of the parallelogram corners, now using datetime values
# # #   parallelogram_data <- tibble(
# # #     x = c(lower_left, lower_right, upper_right, upper_left, lower_left),
# # #     y = c(lower_left[2], lower_right[2], upper_right[2], upper_left[2], lower_left[2])
# # #   )
# # #
# # #   # Create the base plot with the full profile line
# # #   plot <- ggplot2::ggplot(profile_data, ggplot2::aes(x = .data$datetime, y = .data$melatonin)) +
# # #     ggplot2::geom_line(
# # #       ggplot2::aes(color = "Full Profile", group = 1, linetype = "Full Profile"),
# # #       size = 1.25
# # #     ) +
# # #     ggplot2::scale_x_datetime(
# # #       labels = scales::date_format("%H:%M"),
# # #       date_breaks = "2 hours"
# # #     ) +
# # #     ggplot2::labs(
# # #       title = "Melatonin Profile with Parallelogram",
# # #       x = "Time",
# # #       y = "Melatonin Concentration"
# # #     ) +
# # #     ggplot2::theme_minimal() +
# # #     ggplot2::theme(legend.position = "right") +
# # #     ggplot2::scale_linetype_manual(values = c("Full Profile" = "dotted"))
# # #
# # #   # Add base and ascending segments
# # #   plot <- plot +
# # #     ggplot2::geom_line(
# # #       data = dplyr::filter(profile_data, .data$base == 1),
# # #       ggplot2::aes(x = .data$datetime, y = .data$melatonin, color = "Base Segment"),
# # #       size = 1.25
# # #     ) +
# # #     ggplot2::geom_line(
# # #       data = dplyr::filter(profile_data, .data$ascending == 1),
# # #       ggplot2::aes(x = .data$datetime, y = .data$melatonin, color = "Ascending Segment"),
# # #       size = 1.25
# # #     )
# # #
# # #   # Plot the parallelogram
# # #   plot <- plot +
# # #     ggplot2::geom_polygon(
# # #       data = parallelogram_data,
# # #       ggplot2::aes(x = .data$x, y = .data$y),
# # #       fill = "blue", alpha = 0.2, color = "black"
# # #     )
# # #
# # #   # Use annotate to plot the diagonals of the parallelogram
# # #   plot <- plot +
# # #     ggplot2::annotate("segment",
# # #                       x = lower_left[1], y = lower_left[2],
# # #                       xend = upper_right[1], yend = upper_right[2],
# # #                       color = "red", linetype = "dashed", size = 1
# # #     ) +
# # #     ggplot2::annotate("segment",
# # #                       x = lower_right[1], y = lower_right[2],
# # #                       xend = upper_left[1], yend = upper_left[2],
# # #                       color = "red", linetype = "dashed", size = 1
# # #     )
# # #
# # #   return(plot)
# # # }
# #
# # # FROM HERE ON
# # # # Plot the profile and the parallelogram
# # # plot_parallelogram <- function(profile_data, params) {
# # #   # Step 1: Convert datetime to numeric (seconds since Unix epoch)
# # #   profile_data$datetime_numeric <- as.numeric(profile_data$datetime)
# # #
# # #   # Step 2: Perform optimization using optimize_parallelogram
# # #   # Assuming params contains the starting values for the optimization
# # #   opt_result <- optimize_parallelogram(profile_data, params)
# # #
# # #   # The result will be the optimized parameters, e.g.:
# # #   # opt_result[1] = param1, opt_result[2] = param2, opt_result[3] = param3
# # #
# # #   # Step 3: Get the corners of the parallelogram with the optimized parameters
# # #   y0 <- min(profile_data$melatonin)
# # #   y1 <- max(profile_data$melatonin)
# # #   corners <- get_corners(opt_result[1], y0, opt_result[2], y1, opt_result[3])
# # #
# # #   # Extract the corners (they are in numeric form at this point)
# # #   lower_left <- corners[[1]]
# # #   lower_right <- corners[[2]]
# # #   upper_right <- corners[[3]]
# # #   upper_left <- corners[[4]]
# # #
# # #   # Step 4: Convert the numeric values of corners back to datetime
# # #   start_datetime <- min(profile_data$datetime)  # Reference datetime (start)
# # #
# # #   # Convert numeric seconds back to datetime
# # #   lower_left_datetime <- start_datetime + as.difftime(lower_left[1], units = "secs")
# # #   lower_right_datetime <- start_datetime + as.difftime(lower_right[1], units = "secs")
# # #   upper_right_datetime <- start_datetime + as.difftime(upper_right[1], units = "secs")
# # #   upper_left_datetime <- start_datetime + as.difftime(upper_left[1], units = "secs")
# # #
# # #   # Step 5: Create the parallelogram data for plotting
# # #   parallelogram_data <- data.frame(
# # #     x = c(lower_left_datetime, lower_right_datetime, upper_right_datetime, upper_left_datetime),
# # #     y = c(lower_left[2], lower_right[2], upper_right[2], upper_left[2])
# # #   )
# # #
# # #   # Step 6: Plot the profile and overlay the parallelogram
# # #   plot <- ggplot2::ggplot(profile_data, ggplot2::aes(x = .data$datetime, y = .data$melatonin)) +
# # #     # Plot the full profile
# # #     ggplot2::geom_line(
# # #       ggplot2::aes(color = "Full Profile", group = 1, linetype = "Full Profile"),
# # #       size = 1.25
# # #     ) +
# # #     # Overlay the parallelogram
# # #     ggplot2::geom_polygon(
# # #       data = parallelogram_data,
# # #       ggplot2::aes(x = x, y = y),
# # #       fill = "red",
# # #       alpha = 0.5
# # #     ) +
# # #     # Format the x-axis to show only time
# # #     ggplot2::scale_x_datetime(
# # #       labels = scales::date_format("%H:%M"),
# # #       date_breaks = "2 hours"
# # #     ) +
# # #     # Add plot labels
# # #     ggplot2::labs(
# # #       title = "Melatonin Profile with Parallelogram Overlay",
# # #       x = "Time",
# # #       y = "Melatonin Concentration"
# # #     ) +
# # #     ggplot2::theme_minimal() +
# # #     ggplot2::theme(legend.position = "right") +
# # #     ggplot2::scale_linetype_manual(values = c("Full Profile" = "dotted"))
# # #
# # #   # Return the plot
# # #   return(plot)
# # # }
# # #
# # #
# # #
# # # # Helper function to compute parallelogram corners
# # # get_corners <- function(x0, y0, x1, y1, slope) {
# # #   # Ensure x0 and x1 are numeric (in seconds since epoch)
# # #   x0_numeric <- as.numeric(x0)
# # #   x1_numeric <- as.numeric(x1)
# # #
# # #   # Calculate the parallelogram corners using numeric time values
# # #   if (slope == 0) {
# # #     lower_left <- c(x0_numeric, y0)
# # #     lower_right <- c(x1_numeric, y0)
# # #     upper_left <- c(x0_numeric, y1)
# # #     upper_right <- c(x1_numeric, y1)
# # #   } else {
# # #     height <- y1 - y0
# # #     delta_x <- height / slope
# # #
# # #     if (abs(slope) > 1e3) {
# # #       delta_x <- 0
# # #     }
# # #
# # #     lower_left <- c(x0_numeric, y0)
# # #     lower_right <- c(x1_numeric, y0)
# # #     upper_left <- c(x0_numeric + delta_x, y1)
# # #     upper_right <- c(x1_numeric + delta_x, y1)
# # #   }
# # #
# # #   # Return the corners as numeric values (in seconds since epoch)
# # #   return(list(lower_left, lower_right, upper_right, upper_left))
# # # }
# # #
# # #
# # #
# # #
# # #
# # #
# # # # Helper function to compute constraints
# # # constraints <- function(params, x, y, y0, y1) {
# # #   x0 <- params[1]
# # #   x1 <- params[2]
# # #   slope <- params[3]
# # #   constraint_vals <- numeric()
# # #
# # #   corners <- get_corners(x0, y0, x1, y1, slope)
# # #   lower_left <- corners[[1]]
# # #   lower_right <- corners[[2]]
# # #   upper_right <- corners[[3]]
# # #   upper_left <- corners[[4]]
# # #
# # #   for (i in seq_along(x)) {
# # #     xi <- x[i]
# # #     yi <- y[i]
# # #
# # #     y_lower <- if (xi <= lower_right[1]) y0 else y0 + slope * (xi - lower_right[1])
# # #     y_upper <- if (xi >= upper_left[1]) y1 else y1 + slope * (xi - upper_left[1])
# # #
# # #     x_lower <- min(lower_left[1], upper_left[1])
# # #     x_upper <- max(upper_right[1], lower_right[1])
# # #
# # #     constraint_vals <- c(constraint_vals, y_upper - yi)
# # #     constraint_vals <- c(constraint_vals, yi - y_lower)
# # #     constraint_vals <- c(constraint_vals, xi - x_lower)
# # #     constraint_vals <- c(constraint_vals, x_upper - xi)
# # #   }
# # #   return(constraint_vals)
# # # }
# # #
# # # # Helper function for objective minimization
# # # objective <- function(params, x, y, y0, y1) {
# # #   x0 <- params[1]
# # #   x1 <- params[2]
# # #   slope <- params[3]
# # #   corners <- get_corners(x0, y0, x1, y1, slope)
# # #   v1 <- corners[[2]] - corners[[1]]
# # #   v2 <- corners[[4]] - corners[[1]]
# # #   area <- abs(v1[1] * v2[2] - v1[2] * v2[1])
# # #
# # #   c_penalty <- min(constraints(params, x, y, y0, y1))
# # #   if (c_penalty > 0) {
# # #     c_penalty <- 0
# # #   }
# # #   return(area + 1e3 * c_penalty^2)
# # # }
# # #
# # # # Function to optimize parallelogram using indices
# # # optimize_parallelogram <- function(profile_data) {
# # #   # Extract ascending segment and convert time to indices
# # #   ascending_data <- profile_data %>% dplyr::filter(.data$ascending == 1)
# # #
# # #   if (nrow(ascending_data) < 2) {
# # #     warning("Insufficient ascending segment to apply parallelogram rule.")
# # #     return(NULL)
# # #   }
# # #
# # #   x <- seq_along(ascending_data$melatonin)  # Use indices instead of datetime values
# # #   y <- ascending_data$melatonin
# # #
# # #   # Initial guess for x0, x1, and slope
# # #   y0 <- min(y)
# # #   y1 <- max(y)
# # #   x0_initial <- 1
# # #   x1_initial <- length(x)
# # #   slope_initial <- (y1 - y0) / (x1_initial - x0_initial)
# # #
# # #   initial_guess <- c(x0_initial, x1_initial, slope_initial)
# # #
# # #   result <- optim(
# # #     par = initial_guess,
# # #     fn = objective,
# # #     method = "L-BFGS-B",
# # #     lower = c(1, 1, -Inf),
# # #     upper = c(length(x), length(x), Inf),
# # #     x = x, y = y, y0 = y0, y1 = y1
# # #   )
# # #
# # #   return(result$par)
# # # }
# # #
# # # # Function to calculate diagonal slopes of the parallelogram
# # # calculate_diagonal_slopes <- function(params, y0, y1) {
# # #   x0 <- params[1]
# # #   x1 <- params[2]
# # #   slope <- params[3]
# # #
# # #   corners <- get_corners(x0, y0, x1, y1, slope)
# # #   lower_left <- corners[[1]]
# # #   upper_right <- corners[[3]]
# # #   upper_left <- corners[[4]]
# # #   lower_right <- corners[[2]]
# # #
# # #   # Long diagonal
# # #   long_diag_slope <- (upper_right[2] - lower_left[2]) / (upper_right[1] - lower_left[1])
# # #
# # #   # Short diagonal
# # #   short_diag_slope <- (lower_right[2] - upper_left[2]) / (lower_right[1] - upper_left[1])
# # #
# # #   return(c(long_diag_slope = long_diag_slope, short_diag_slope = short_diag_slope))
# # # }
# # #
# # # # Main function to truncate by parallelogram rule
# # # truncate_by_parallelogram_rule <- function(profile_data) {
# # #   # Extract ascending segment
# # #   ascending_data <- profile_data %>% dplyr::filter(.data$ascending == 1)
# # #
# # #   if (nrow(ascending_data) < 2) {
# # #     warning("Insufficient ascending segment to apply parallelogram rule.")
# # #     return(profile_data)
# # #   }
# # #
# # #   # Optimize the parallelogram
# # #   params <- optimize_parallelogram(profile_data)
# # #   print(params)
# # #
# # #   if (is.null(params)) return(profile_data)
# # #
# # #   # Calculate diagonal slopes
# # #   y0 <- min(ascending_data$melatonin)
# # #   y1 <- max(ascending_data$melatonin)
# # #   diagonals <- calculate_diagonal_slopes(params, y0, y1)
# # #   long_diag_slope <- diagonals["long_diag_slope"]
# # #   short_diag_slope <- diagonals["short_diag_slope"]
# # #
# # #   # Check the slope condition
# # #   if (abs(long_diag_slope) < 0.5 * abs(short_diag_slope)) {
# # #     # Truncate the ascending segment
# # #     profile_data <- profile_data %>%
# # #       dplyr::mutate(
# # #         ascending = ifelse(.data$ascending == 1 & .data$datetime > ascending_data$datetime[1], 0, .data$ascending)
# # #       )
# # #   }
# # #
# # #   # Plot the profile with parallelogram and diagonals
# # #   plot <- plot_parallelogram(profile_data, params)
# # #   print(plot)
# # #   return(profile_data)
# # # }
# # #
# #
# # # # Function to perform the parallelogram optimization
# # # optimize_parallelogram <- function(profile_data) {
# # #   # Extract ascending segment and convert time to indices
# # #   ascending_data <- profile_data %>% dplyr::filter(.data$ascending == 1)
# # #
# # #   if (nrow(ascending_data) < 2) {
# # #     warning("Insufficient ascending segment to apply parallelogram rule.")
# # #     return(NULL)
# # #   }
# # #
# # #   x <- seq_along(ascending_data$melatonin)  # Use indices instead of datetime values
# # #   y <- ascending_data$melatonin
# # #
# # #   # Initial guess for x0, x1, and slope
# # #   y0 <- min(y)
# # #   y1 <- max(y)
# # #   x0_initial <- 1
# # #   x1_initial <- length(x)
# # #   slope_initial <- (y1 - y0) / (x1_initial - x0_initial)
# # #
# # #   initial_guess <- c(x0_initial, x1_initial, slope_initial)
# # #
# # #   result <- optim(
# # #     par = initial_guess,
# # #     fn = objective,
# # #     method = "L-BFGS-B",
# # #     lower = c(1, 1, -Inf),
# # #     upper = c(length(x), length(x), Inf),
# # #     x = x, y = y, y0 = y0, y1 = y1
# # #   )
# # #
# # #   return(result$par)  # Return optimized parameters (x0, x1, slope)
# # # }
# # #
# # # # Function to calculate diagonal slopes of the parallelogram
# # # calculate_diagonal_slopes <- function(params, y0, y1) {
# # #   x0 <- params[1]
# # #   x1 <- params[2]
# # #   slope <- params[3]
# # #
# # #   corners <- get_corners(x0, y0, x1, y1, slope)
# # #   lower_left <- corners[[1]]
# # #   upper_right <- corners[[3]]
# # #   upper_left <- corners[[4]]
# # #   lower_right <- corners[[2]]
# # #
# # #   # Long diagonal
# # #   long_diag_slope <- (upper_right[2] - lower_left[2]) / (upper_right[1] - lower_left[1])
# # #
# # #   # Short diagonal
# # #   short_diag_slope <- (lower_right[2] - upper_left[2]) / (lower_right[1] - upper_left[1])
# # #
# # #   return(c(long_diag_slope = long_diag_slope, short_diag_slope = short_diag_slope))
# # # }
# # #
# # # # Function to calculate the parallelogram corners
# # # get_corners <- function(x0, y0, x1, y1, slope) {
# # #   # Ensure x0 and x1 are numeric (in seconds since epoch)
# # #   x0_numeric <- as.numeric(x0)
# # #   x1_numeric <- as.numeric(x1)
# # #
# # #   # Calculate the parallelogram corners using numeric time values
# # #   if (slope == 0) {
# # #     lower_left <- c(x0_numeric, y0)
# # #     lower_right <- c(x1_numeric, y0)
# # #     upper_left <- c(x0_numeric, y1)
# # #     upper_right <- c(x1_numeric, y1)
# # #   } else {
# # #     height <- y1 - y0
# # #     delta_x <- height / slope
# # #
# # #     if (abs(slope) > 1e3) {
# # #       delta_x <- 0
# # #     }
# # #
# # #     lower_left <- c(x0_numeric, y0)
# # #     lower_right <- c(x1_numeric, y0)
# # #     upper_left <- c(x0_numeric + delta_x, y1)
# # #     upper_right <- c(x1_numeric + delta_x, y1)
# # #   }
# # #
# # #   # Return the corners as numeric values (in seconds since epoch)
# # #   return(list(lower_left, lower_right, upper_right, upper_left))
# # # }
# # #
# # # # Function to compute constraints (used in the optimization)
# # # constraints <- function(params, x, y, y0, y1) {
# # #   x0 <- params[1]
# # #   x1 <- params[2]
# # #   slope <- params[3]
# # #   constraint_vals <- numeric()
# # #
# # #   corners <- get_corners(x0, y0, x1, y1, slope)
# # #   lower_left <- corners[[1]]
# # #   lower_right <- corners[[2]]
# # #   upper_right <- corners[[3]]
# # #   upper_left <- corners[[4]]
# # #
# # #   for (i in seq_along(x)) {
# # #     xi <- x[i]
# # #     yi <- y[i]
# # #
# # #     y_lower <- if (xi <= lower_right[1]) y0 else y0 + slope * (xi - lower_right[1])
# # #     y_upper <- if (xi >= upper_left[1]) y1 else y1 + slope * (xi - upper_left[1])
# # #
# # #     x_lower <- min(lower_left[1], upper_left[1])
# # #     x_upper <- max(upper_right[1], lower_right[1])
# # #
# # #     constraint_vals <- c(constraint_vals, y_upper - yi)
# # #     constraint_vals <- c(constraint_vals, yi - y_lower)
# # #     constraint_vals <- c(constraint_vals, xi - x_lower)
# # #     constraint_vals <- c(constraint_vals, x_upper - xi)
# # #   }
# # #   return(constraint_vals)
# # # }
# # #
# # # # Function to compute the objective function for optimization
# # # objective <- function(params, x, y, y0, y1) {
# # #   x0 <- params[1]
# # #   x1 <- params[2]
# # #   slope <- params[3]
# # #   corners <- get_corners(x0, y0, x1, y1, slope)
# # #   v1 <- corners[[2]] - corners[[1]]
# # #   v2 <- corners[[4]] - corners[[1]]
# # #   area <- abs(v1[1] * v2[2] - v1[2] * v2[1])
# # #
# # #   c_penalty <- min(constraints(params, x, y, y0, y1))
# # #   if (c_penalty > 0) {
# # #     c_penalty <- 0
# # #   }
# # #   return(area + 1e3 * c_penalty^2)
# # # }
# # #
# # # # Main function to truncate by parallelogram rule
# # # truncate_by_parallelogram_rule <- function(profile_data) {
# # #   # Extract ascending segment
# # #   ascending_data <- profile_data %>% dplyr::filter(.data$ascending == 1)
# # #
# # #   if (nrow(ascending_data) < 2) {
# # #     warning("Insufficient ascending segment to apply parallelogram rule.")
# # #     return(profile_data)
# # #   }
# # #
# # #   # Optimize the parallelogram
# # #   params <- optimize_parallelogram(profile_data)
# # #   print(params)
# # #
# # #   if (is.null(params)) return(profile_data)
# # #
# # #   # Calculate diagonal slopes
# # #   y0 <- min(ascending_data$melatonin)
# # #   y1 <- max(ascending_data$melatonin)
# # #   diagonals <- calculate_diagonal_slopes(params, y0, y1)
# # #   long_diag_slope <- diagonals["long_diag_slope"]
# # #   short_diag_slope <- diagonals["short_diag_slope"]
# # #
# # #   # Check the slope condition
# # #   if (abs(long_diag_slope) < 0.5 * abs(short_diag_slope)) {
# # #     # Truncate the ascending segment
# # #     profile_data <- profile_data %>%
# # #       dplyr::mutate(
# # #         ascending = ifelse(.data$ascending == 1 & .data$datetime > ascending_data$datetime[1], 0, .data$ascending)
# # #       )
# # #   }
# # #
# # #   # Plot the profile with parallelogram and diagonals
# # #   plot <- plot_parallelogram(profile_data, params)  # Pass optimized params to the plot
# # #   print(plot)
# # #
# # #   return(profile_data)
# # # }
# # #
# # # # Plotting function to overlay the parallelogram
# # # plot_parallelogram <- function(profile_data, params) {
# # #   # Step 1: Convert datetime to numeric (seconds since Unix epoch)
# # #   profile_data$datetime_numeric <- as.numeric(profile_data$datetime)
# # #
# # #   # Step 2: Get the corners of the parallelogram with the optimized parameters
# # #   y0 <- min(profile_data$melatonin)
# # #   y1 <- max(profile_data$melatonin)
# # #   corners <- get_corners(params[1], y0, params[2], y1, params[3])
# # #
# # #   # Extract the corners (they are in numeric form at this point)
# # #   lower_left <- corners[[1]]
# # #   lower_right <- corners[[2]]
# # #   upper_right <- corners[[3]]
# # #   upper_left <- corners[[4]]
# # #
# # #   # Step 3: Convert the numeric values of corners back to datetime
# # #   start_datetime <- min(profile_data$datetime)  # Reference datetime (start)
# # #   lower_left_datetime <- start_datetime + as.difftime(lower_left[1], units = "secs")
# # #   lower_right_datetime <- start_datetime + as.difftime(lower_right[1], units = "secs")
# # #   upper_right_datetime <- start_datetime + as.difftime(upper_right[1], units = "secs")
# # #   upper_left_datetime <- start_datetime + as.difftime(upper_left[1], units = "secs")
# # #
# # #   # Step 4: Create the parallelogram data for plotting
# # #   parallelogram_data <- data.frame(
# # #     x = c(lower_left_datetime, lower_right_datetime, upper_right_datetime, upper_left_datetime),
# # #     y = c(lower_left[2], lower_right[2], upper_right[2], upper_left[2])
# # #   )
# # #
# # #   # Step 5: Plot the profile and overlay the parallelogram
# # #   plot <- ggplot2::ggplot(profile_data, ggplot2::aes(x = .data$datetime, y = .data$melatonin)) +
# # #     # Plot the full profile
# # #     ggplot2::geom_line(
# # #       ggplot2::aes(color = "Full Profile", group = 1, linetype = "Full Profile"),
# # #       size = 1.25
# # #     ) +
# # #     # Overlay the parallelogram
# # #     ggplot2::geom_polygon(
# # #       data = parallelogram_data,
# # #       ggplot2::aes(x = x, y = y),
# # #       fill = "red",
# # #       alpha = 0.5
# # #     ) +
# # #     # Format the x-axis to show only time
# # #     ggplot2::scale_x_datetime(
# # #       labels = scales::date_format("%H:%M"),
# # #       date_breaks = "2 hours"
# # #     ) +
# # #     # Add plot labels
# # #     ggplot2::labs(
# # #       title = "Melatonin Profile with Parallelogram Overlay",
# # #       x = "Time",
# # #       y = "Melatonin Concentration"
# # #     ) +
# # #     ggplot2::theme_minimal() +
# # #     ggplot2::theme(legend.position = "right") +
# # #     ggplot2::scale_linetype_manual(values = c("Full Profile" = "dotted"))
# # #
# # #   # Return the plot
# # #   return(plot)
# # # }
# # #
# #
# #
# # # Truncation and Optimization
# # ########################################
# # # Function to optimize parallelogram
# # optimize_parallelogram <- function(profile_data) {
# #   # Extract ascending segment
# #   ascending_data <- profile_data %>% dplyr::filter(.data$ascending == 1)
# #
# #   if (nrow(ascending_data) < 2) {
# #     warning("Insufficient ascending segment to apply parallelogram rule.")
# #     return(NULL)
# #   }
# #
# #   # Use indices instead of datetime for optimization
# #   x <- seq_along(ascending_data$melatonin)  # Indices of ascending data
# #   y <- ascending_data$melatonin
# #
# #   # Initial guesses
# #   y0 <- min(y)
# #   y1 <- max(y)
# #   x0_initial <- 1
# #   x1_initial <- length(x)
# #   slope_initial <- (y1 - y0) / (x1_initial - x0_initial)
# #
# #   initial_guess <- c(x0_initial, x1_initial, slope_initial)
# #
# #   # Perform optimization
# #   result <- optim(
# #     par = initial_guess,
# #     fn = objective,
# #     method = "L-BFGS-B",
# #     lower = c(1, 1, -Inf),
# #     upper = c(length(x), length(x), Inf),
# #     x = x, y = y, y0 = y0, y1 = y1
# #   )
# #
# #   return(result$par)  # Optimized parameters
# # }
# #
# # # Helper function for objective minimization
# # objective <- function(params, x, y, y0, y1) {
# #   x0 <- params[1]
# #   x1 <- params[2]
# #   slope <- params[3]
# #   corners <- get_corners(x0, y0, x1, y1, slope)
# #
# #   # Extract corners
# #   v1 <- corners[[2]] - corners[[1]]
# #   v2 <- corners[[4]] - corners[[1]]
# #
# #   # Calculate the area of the parallelogram
# #   area <- abs(v1[1] * v2[2] - v1[2] * v2[1])
# #
# #   # Constraints penalty
# #   c_penalty <- min(constraints(params, x, y, y0, y1))
# #   if (c_penalty > 0) {
# #     c_penalty <- 0
# #   }
# #
# #   # Return the objective value (area + penalty for violations)
# #   return(area + 1e3 * c_penalty^2)
# # }
# #
# # constraints <- function(params, x, y, y0, y1) {
# #   x0 <- params[1]
# #   x1 <- params[2]
# #   slope <- params[3]
# #   constraint_vals <- numeric()
# #
# #   corners <- get_corners(x0, y0, x1, y1, slope)
# #   lower_left <- corners[[1]]
# #   lower_right <- corners[[2]]
# #   upper_right <- corners[[3]]
# #   upper_left <- corners[[4]]
# #
# #   for (i in seq_along(x)) {
# #     xi <- x[i]
# #     yi <- y[i]
# #
# #     y_lower <- if (xi <= lower_right[1]) y0 else y0 + slope * (xi - lower_right[1])
# #     y_upper <- if (xi >= upper_left[1]) y1 else y1 + slope * (xi - upper_left[1])
# #
# #     x_lower <- min(lower_left[1], upper_left[1])
# #     x_upper <- max(upper_right[1], lower_right[1])
# #
# #     constraint_vals <- c(constraint_vals, y_upper - yi)
# #     constraint_vals <- c(constraint_vals, yi - y_lower)
# #     constraint_vals <- c(constraint_vals, xi - x_lower)
# #     constraint_vals <- c(constraint_vals, x_upper - xi)
# #   }
# #   return(constraint_vals)
# # }
# #
# # # Function to calculate diagonal slopes of the parallelogram
# # calculate_diagonal_slopes <- function(params, y0, y1) {
# #   x0 <- params[1]
# #   x1 <- params[2]
# #   slope <- params[3]
# #
# #   # Get the parallelogram corners
# #   corners <- get_corners(x0, y0, x1, y1, slope)
# #
# #   # Extract corners
# #   lower_left <- corners[[1]]
# #   lower_right <- corners[[2]]
# #   upper_right <- corners[[3]]
# #   upper_left <- corners[[4]]
# #
# #   # Calculate the slopes of the diagonals
# #   # Long diagonal (from lower-left to upper-right)
# #   long_diag_slope <- (upper_right[2] - lower_left[2]) / (upper_right[1] - lower_left[1])
# #
# #   # Short diagonal (from upper-left to lower-right)
# #   short_diag_slope <- (lower_right[2] - upper_left[2]) / (lower_right[1] - upper_left[1])
# #
# #   # Return both slopes
# #   return(c(long_diag_slope = long_diag_slope, short_diag_slope = short_diag_slope))
# # }
# #
# # get_corners <- function(x0, y0, x1, y1, slope) {
# #   if (slope == 0) {
# #     lower_left <- c(x0, y0)
# #     lower_right <- c(x1, y0)
# #     upper_left <- c(x0, y1)
# #     upper_right <- c(x1, y1)
# #   } else {
# #     height <- y1 - y0
# #     delta_x <- height / slope
# #     if (abs(slope) > 1e3) {
# #       delta_x <- 0
# #     }
# #
# #     lower_left <- c(x0, y0)
# #     lower_right <- c(x1, y0)
# #     upper_left <- c(x0 + delta_x, y1)
# #     upper_right <- c(x1 + delta_x, y1)
# #   }
# #   return(list(lower_left, lower_right, upper_right, upper_left))
# # }
# #
# # # Truncation logic
# # truncate_by_parallelogram_rule <- function(profile_data) {
# #   # Extract ascending segment
# #   ascending_data <- profile_data %>% dplyr::filter(.data$ascending == 1)
# #
# #   if (nrow(ascending_data) < 2) {
# #     warning("Insufficient ascending segment to apply parallelogram rule.")
# #     return(profile_data)
# #   }
# #
# #   # Optimize the parallelogram
# #   params <- optimize_parallelogram(profile_data)
# #   if (is.null(params)) {
# #     warning("Optimization failed.")
# #     return(profile_data)
# #   }
# #
# #   # Print optimized parameters
# #   cat("Optimized Parameters:\n")
# #   print(params)
# #
# #   # Plot the profile with the parallelogram overlay
# #   plot <- plot_profile_with_regions_and_parallelogram(profile_data, params)
# #   print(plot)
# #
# #   return(profile_data)
# # }
# #
# #
# #
# # ########################################
# # # Plotting the Profile and Parallelogram
# # plot_profile_with_regions_and_parallelogram <- function(profile_data, params) {
# #   # Step 1: Calculate the parallelogram corners
# #   y0 <- min(profile_data$melatonin)
# #   y1 <- max(profile_data$melatonin)
# #
# #   # Debug: Print Y-axis bounds
# #   cat("Y-Axis Bounds: y0 =", y0, ", y1 =", y1, "\n")
# #
# #   corners <- get_corners(params[1], y0, params[2], y1, params[3])
# #
# #   lower_left <- corners[[1]]
# #   lower_right <- corners[[2]]
# #   upper_right <- corners[[3]]
# #   upper_left <- corners[[4]]
# #
# #   # Debug: Print corner coordinates (numeric)
# #   cat("Parallelogram Corners (Numeric):\n")
# #   print(list(lower_left, lower_right, upper_right, upper_left))
# #
# #   # Step 2: Convert numeric corners to datetime for plotting
# #   start_datetime <- min(profile_data$datetime)
# #   parallelogram_data <- data.frame(
# #     x = c(
# #       start_datetime + as.difftime(lower_left[1], units = "secs"),
# #       start_datetime + as.difftime(lower_right[1], units = "secs"),
# #       start_datetime + as.difftime(upper_right[1], units = "secs"),
# #       start_datetime + as.difftime(upper_left[1], units = "secs")
# #     ),
# #     y = c(lower_left[2], lower_right[2], upper_right[2], upper_left[2])
# #   )
# #
# #   # Debug: Print datetime-converted parallelogram corners
# #   cat("Parallelogram Corners (Datetime):\n")
# #   print(parallelogram_data)
# #
# #
# #   # Step 3: Plot the profile with base, ascending regions, and parallelogram overlay
# #   plot <- ggplot2::ggplot(profile_data, ggplot2::aes(x = .data$datetime, y = .data$melatonin)) +
# #     ggplot2::geom_line(ggplot2::aes(color = "Full Profile", group = 1), size = 1.25) +
# #     ggplot2::geom_line(
# #       data = dplyr::filter(profile_data, base == 1),
# #       ggplot2::aes(color = "Base Segment"),
# #       size = 1.25
# #     ) +
# #     ggplot2::geom_line(
# #       data = dplyr::filter(profile_data, ascending == 1),
# #       ggplot2::aes(color = "Ascending Segment"),
# #       size = 1.25
# #     ) +
# #     ggplot2::geom_polygon(
# #       data = parallelogram_data,
# #       ggplot2::aes(x = x, y = y, group = 1),
# #       fill = "red",
# #       alpha = 0.5
# #     ) +
# #     ggplot2::scale_x_datetime(labels = scales::date_format("%H:%M"), date_breaks = "2 hours") +
# #     ggplot2::labs(title = "Melatonin Profile with Parallelogram Overlay",
# #                   x = "Time", y = "Melatonin Concentration") +
# #     ggplot2::theme_minimal() +
# #     ggplot2::scale_color_manual(
# #       values = c("Full Profile" = "black", "Base Segment" = "blue", "Ascending Segment" = "green"),
# #       name = "Segments"
# #     ) +
# #     ggplot2::theme(legend.position = "right")
# #
# #   return(plot)
# #
# # }
# #
# # get_corners <- function(x0, y0, x1, y1, slope) {
# #   # Ensure x0 and x1 are numeric
# #   x0_numeric <- as.numeric(x0)
# #   x1_numeric <- as.numeric(x1)
# #
# #   # Debug: Print input parameters for corner calculation
# #   cat("Inputs to get_corners:\n")
# #   cat("x0 =", x0_numeric, "x1 =", x1_numeric, "y0 =", y0, "y1 =", y1, "slope =", slope, "\n")
# #
# #   # Calculate corners
# #   if (slope == 0) {
# #     lower_left <- c(x0_numeric, y0)
# #     lower_right <- c(x1_numeric, y0)
# #     upper_left <- c(x0_numeric, y1)
# #     upper_right <- c(x1_numeric, y1)
# #   } else {
# #     height <- y1 - y0
# #     delta_x <- height / slope
# #
# #     lower_left <- c(x0_numeric, y0)
# #     lower_right <- c(x1_numeric, y0)
# #     upper_left <- c(x0_numeric + delta_x, y1)
# #     upper_right <- c(x1_numeric + delta_x, y1)
# #   }
# #
# #   # Debug: Print calculated corners
# #   cat("Calculated Corners (Numeric):\n")
# #   print(data.frame(
# #     x = c(lower_left[1], lower_right[1], upper_right[1], upper_left[1]),
# #     y = c(lower_left[2], lower_right[2], upper_right[2], upper_left[2])
# #   ))
# #
# #   return(list(lower_left, lower_right, upper_right, upper_left))
# # }
# #
# # optimize_parallelogram <- function(profile_data) {
# #   # Extract ascending segment
# #   ascending_data <- profile_data %>% dplyr::filter(.data$ascending == 1)
# #
# #   if (nrow(ascending_data) < 2) {
# #     warning("Insufficient ascending segment to optimize parallelogram.")
# #     return(NULL)
# #   }
# #
# #   x <- seq_along(ascending_data$melatonin)  # Use indices instead of datetime values
# #   y <- ascending_data$melatonin
# #
# #   # Initial guess for x0, x1, and slope
# #   y0 <- min(y)
# #   y1 <- max(y)
# #   x0_initial <- 1
# #   x1_initial <- length(x)
# #   slope_initial <- (y1 - y0) / (x1_initial - x0_initial)
# #
# #   initial_guess <- c(x0_initial, x1_initial, slope_initial)
# #
# #   result <- optim(
# #     par = initial_guess,
# #     fn = objective,
# #     method = "L-BFGS-B",
# #     lower = c(1, 1, -Inf),
# #     upper = c(length(x), length(x), Inf),
# #     x = x, y = y, y0 = y0, y1 = y1
# #   )
# #
# #   # Debug: Print optimization result
# #   cat("Optimization Result:\n")
# #   print(result$par)
# #
# #   return(result$par)
# # }
# #
#
# library(ggplot2)
# library(dplyr)
#
# # Helper function to convert POSIXct to numeric (seconds since the epoch)
# datetime_to_numeric <- function(datetime) {
#   as.numeric(datetime)
# }
#
# # Function to compute the parallelogram corners based on the optimized parameters
# get_corners <- function(x1, x2, y0, y1, slope) {
#   # Calculate the parallelogram corners using numeric time values
#   delta_x <- (y1 - y0) / slope  # Calculate horizontal shift based on slope
#   lower_left <- c(x1, y0)
#   lower_right <- c(x2, y0)
#   upper_left <- c(x1 + delta_x, y1)
#   upper_right <- c(x2 + delta_x, y1)
#
#   return(list(lower_left, lower_right, upper_right, upper_left))
# }
#
# # Objective function for optimization (simplified)
# objective <- function(params, x, y, y0, y1) {
#   # Parameters: x1, x2, slope
#   x1 <- params[1]
#   x2 <- params[2]
#   slope <- params[3]
#
#   # Calculate the corners of the parallelogram
#   corners <- get_corners(x1, x2, y0, y1, slope)
#
#   # Calculate the area of the parallelogram using cross-product of vectors
#   v1 <- corners[[2]] - corners[[1]]
#   v2 <- corners[[4]] - corners[[1]]
#   area <- abs(v1[1] * v2[2] - v1[2] * v2[1])
#
#   return(area)
# }
#
# # Simplified Optimization function
# optimize_parallelogram <- function(profile_data) {
#   # Extract ascending segment and convert time to numeric
#   ascending_data <- profile_data %>% filter(ascending == 1)
#   if (nrow(ascending_data) < 2) {
#     warning("Insufficient ascending segment to apply parallelogram rule.")
#     return(NULL)
#   }
#
#   # Convert datetime to numeric (seconds since Unix epoch)
#   x <- datetime_to_numeric(ascending_data$datetime)
#   y <- ascending_data$melatonin
#
#   # Initial guess for x1, x2, and slope
#   y0 <- min(y)
#   y1 <- max(y)
#   x1_initial <- min(x)
#   x2_initial <- max(x)
#   slope_initial <- (y1 - y0) / (x2_initial - x1_initial)
#   initial_guess <- c(x1_initial, x2_initial, slope_initial)
#
#   # Run optimization to minimize the area (parallelogram area)
#   result <- optim(par = initial_guess, fn = objective, method = "L-BFGS-B",
#                   lower = c(min(x), min(x), -Inf), upper = c(max(x), max(x), Inf),
#                   x = x, y = y, y0 = y0, y1 = y1)
#
#   return(result$par)  # Optimized parameters (x1, x2, slope)
# }
#
# # Function to plot the profile and parallelogram
# plot_parallelogram <- function(profile_data, params) {
#   # Step 1: Convert datetime to numeric (seconds since Unix epoch)
#   profile_data$datetime_numeric <- datetime_to_numeric(profile_data$datetime)
#
#   # Step 2: Get the corners of the parallelogram with the optimized parameters
#   y0 <- min(profile_data$melatonin)
#   y1 <- max(profile_data$melatonin)
#   corners <- get_corners(params[1], params[2], y0, y1, params[3])
#
#   # Extract the corners (they are in numeric form at this point)
#   lower_left <- corners[[1]]
#   lower_right <- corners[[2]]
#   upper_right <- corners[[3]]
#   upper_left <- corners[[4]]
#
#   # Step 3: Convert the numeric values of corners back to datetime
#   start_datetime <- min(profile_data$datetime)  # Reference datetime (start)
#   lower_left_datetime <- as.POSIXct(lower_left[1], origin = "1970-01-01", tz = "UTC")
#   lower_right_datetime <- as.POSIXct(lower_right[1], origin = "1970-01-01", tz = "UTC")
#   upper_right_datetime <- as.POSIXct(upper_right[1], origin = "1970-01-01", tz = "UTC")
#   upper_left_datetime <- as.POSIXct(upper_left[1], origin = "1970-01-01", tz = "UTC")
#
#   # Step 4: Create the parallelogram data for plotting
#   parallelogram_data <- data.frame(
#     x = c(lower_left_datetime, lower_right_datetime, upper_right_datetime, upper_left_datetime),
#     y = c(lower_left[2], lower_right[2], upper_right[2], upper_left[2])
#   )
#
#   # Step 5: Plot the profile and overlay the parallelogram
#   plot <- ggplot(profile_data, aes(x = datetime, y = melatonin)) +
#     # Plot the full profile
#     geom_line(aes(color = "Full Profile", group = 1, linetype = "Full Profile"), size = 1.25) +
#     # Overlay the parallelogram
#     geom_polygon(data = parallelogram_data, aes(x = x, y = y), fill = "red", alpha = 0.5) +
#     # Format the x-axis to show only time
#     scale_x_datetime(labels = scales::date_format("%H:%M"), date_breaks = "2 hours") +
#     # Add plot labels
#     labs(title = "Melatonin Profile with Parallelogram Overlay", x = "Time", y = "Melatonin Concentration") +
#     theme_minimal() +
#     theme(legend.position = "right") +
#     scale_linetype_manual(values = c("Full Profile" = "dotted"))
#
#   return(plot)
# }
#
# # Main function to truncate by parallelogram rule and plot the result
# truncate_by_parallelogram_rule <- function(profile_data) {
#   # Extract ascending segment
#   ascending_data <- profile_data %>% dplyr::filter(ascending == 1)
#
#   if (nrow(ascending_data) < 2) {
#     warning("Insufficient ascending segment to apply parallelogram rule.")
#     return(profile_data)
#   }
#
#   # Optimize the parallelogram
#   params <- optimize_parallelogram(profile_data)
#   print("Optimized Parameters:")
#   print(params)
#
#   if (is.null(params)) return(profile_data)
#
#   # Plot the profile with parallelogram
#   plot <- plot_parallelogram(profile_data, params)
#   print(plot)
#
#   return(profile_data)
# }
#
# # Example:
# # profile_data <- your_profile_data  # Assuming `profile_data` is your dataset
# # profile_data <- truncate_by_parallelogram_rule(profile_data)

# # THIS WORKED 27.11.2024
# # # Function to get corners of the parallelogram
# get_corners <- function(x0, y0, x1, y1, slope) {
#   # x0 <- 19.88825
#   # x1 <- 20.316904
#   # y0 <- min(c(1.661, 9.423, 19.846, 30.492)) # y0 = 1.661
#   # y1 <- max(c(1.661, 9.423, 19.846, 30.492)) # y1 = 30.492
#   # slope <- 9.573235
#   if (slope == 0) {
#     lower_left <- c(x0, y0)
#     lower_right <- c(x1, y0)
#     upper_left <- c(x0, y1)
#     upper_right <- c(x1, y1)
#   } else {
#     height <- y1 - y0
#     delta_x <- height / slope
#
#     if (abs(slope) > 1e3) {
#       delta_x <- 0
#     }
#
#     lower_left <- c(x0, y0)
#     lower_right <- c(x1, y0)
#     upper_left <- c(x0 + delta_x, y1)
#     upper_right <- c(x1 + delta_x, y1)
#   }
#
#   return(list(lower_left, lower_right, upper_right, upper_left))
# }
#
#
# # Function for the constraints
# constraints <- function(params, x, y, y0, y1) {
#   x0 <- params[1]
#   x1 <- params[2]
#   slope <- params[3]
#   constraint_vals <- numeric()
#
#   corners <- get_corners(x0, y0, x1, y1, slope)
#   print("Functional Corners:")
#   print(corners)
#   lower_left <- corners[[1]]
#   lower_right <- corners[[2]]
#   upper_right <- corners[[3]]
#   upper_left <- corners[[4]]
#
#   # Loop over x and y values to compute all constraints
#   for (i in seq_along(x)) {
#     xi <- x[i]
#     yi <- y[i]
#
#     y_lower <- if(xi <= lower_right[1]) y0 else y0 + slope * (xi - lower_right[1])
#     y_upper <- if(xi >= upper_left[1]) y1 else y1 + slope * (xi - upper_left[1])
#
#     x_lower <- min(lower_left[1], upper_left[1])
#     x_upper <- max(upper_right[1], lower_right[1])
#
#     # Append all constraints
#     constraint_vals <- c(
#       constraint_vals,
#       y_upper - yi,
#       yi - y_lower,
#       xi - x_lower,
#       x_upper - xi
#     )
#   }
#
#   # Add debug statement here
#   #print("Functional Constraints:")
#   #print(constraint_vals)
#   constraint_vals
# }

# # Function for the objective
# objective <- function(params, x, y, y0, y1) {
#   x0 <- params[1]
#   x1 <- params[2]
#   slope <- params[3]
#
#   corners <- get_corners(x0, y0, x1, y1, slope)
#   v1 <- corners[[2]] - corners[[1]]
#   v2 <- corners[[4]] - corners[[1]]
#   area <- abs(v1[1] * v2[2] - v1[2] * v2[1])
#
#   c_penalty <- min(constraints(params, x, y, y0, y1))
#   if (c_penalty > 0) {
#     c_penalty <- 0
#   }
#
#   return(area + 1e3 * c_penalty^2)
# }
#
# # Function to optimize parallelogram
# optimize_parallelogram <- function(x, y) {
#   y0 <- min(y)
#   y1 <- max(y)
#   x0_initial <- min(x)
#   x1_initial <- max(x)
#   slope_initial <- (y1 - y0) / (x1_initial - x0_initial)
#
#   initial_guess <- c(x0_initial, x1_initial, slope_initial)
#
#   result <- stats::optim(
#     par = initial_guess,
#     fn = function(params) objective(params, x, y, y0, y1),
#     method = "L-BFGS-B",
#     lower = c(-Inf, -Inf, 0),
#     upper = c(Inf, Inf, Inf)
#   )
#
#   result$par
# }

# # Helper function to calculate corners
# get_corners <- function(x0, y0, x1, y1, slope) {
#   if (slope == 0) {
#     lower_left <- c(x0, y0)
#     lower_right <- c(x1, y0)
#     upper_left <- c(x0, y1)
#     upper_right <- c(x1, y1)
#   } else {
#     height <- y1 - y0
#     delta_x <- height / slope
#
#     if (abs(slope) > 1e3) {
#       delta_x <- 0
#     }
#
#     lower_left <- c(x0, y0)
#     lower_right <- c(x1, y0)
#     upper_left <- c(x0 + delta_x, y1)
#     upper_right <- c(x1 + delta_x, y1)
#   }
#
#   list(lower_left, lower_right, upper_right, upper_left)
# }
#
# # Helper function to calculate constraints
# calculate_constraints <- function(params, x, y, y0, y1) {
#   x0 <- params[1]
#   x1 <- params[2]
#   slope <- params[3]
#
#   corners <- get_corners(x0, y0, x1, y1, slope)
#   lower_left <- corners[[1]]
#   lower_right <- corners[[2]]
#   upper_right <- corners[[3]]
#   upper_left <- corners[[4]]
#
#   purrr::map_dbl(seq_along(x), function(i) {
#     xi <- x[i]
#     yi <- y[i]
#
#     y_lower <- if (xi <= lower_right[1]) y0 else y0 + slope * (xi - lower_right[1])
#     y_upper <- if (xi >= upper_left[1]) y1 else y1 + slope * (xi - upper_left[1])
#
#     x_lower <- min(lower_left[1], upper_left[1])
#     x_upper <- max(upper_right[1], lower_right[1])
#
#     c(
#       y_upper - yi,
#       yi - y_lower,
#       xi - x_lower,
#       x_upper - xi
#     )
#   }) %>% unlist()
# }
#
# # Objective function
# objective_function <- function(params, x, y, y0, y1) {
#   x0 <- params[1]
#   x1 <- params[2]
#   slope <- params[3]
#
#   corners <- get_corners(x0, y0, x1, y1, slope)
#   v1 <- corners[[2]] - corners[[1]]
#   v2 <- corners[[4]] - corners[[1]]
#   area <- abs(v1[1] * v2[2] - v1[2] * v2[1])
#
#   constraints <- calculate_constraints(params, x, y, y0, y1)
#   c_penalty <- min(constraints)
#
#   if (c_penalty > 0) {
#     c_penalty <- 0
#   }
#
#   area + 1e3 * c_penalty^2
# }
#
# # Main optimization function
# optimize_parallelogram <- function(x, y) {
#   y0 <- min(y)
#   y1 <- max(y)
#   x0_initial <- min(x)
#   x1_initial <- max(x)
#   slope_initial <- (y1 - y0) / (x1_initial - x0_initial)
#
#   initial_guess <- c(x0_initial, x1_initial, slope_initial)
#
#   optim_result <- stats::optim(
#     par = initial_guess,
#     fn = function(params) objective_function(params, x, y, y0, y1),
#     method = "L-BFGS-B",
#     lower = c(-Inf, -Inf, 0),
#     upper = c(Inf, Inf, Inf)
#   )
#
#   optim_result$par
# }
