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


plot_parallelogram <- function(plot, profile_data, pll_result) {
  if (is.null(pll_result)) {
    stop("pll_result must be provided to plot the parallelogram.")
  }

  # Extract optimized parameters
  x0_posix <- pll_result$pll_datetime_0
  x1_posix <- pll_result$pll_datetime_1
  slope <- pll_result$pll_slope

  # Filter for ascending data
  profile_data_ascending <- profile_data %>% dplyr::filter(.data$ascending == 1)
  y0 <- min(profile_data_ascending$melatonin)
  y1 <- max(profile_data_ascending$melatonin)

  # Convert datetime to numeric for parallelogram calculations
  x0_numeric <- posixct_to_decimal(x0_posix, profile_data$datetime)
  x1_numeric <- posixct_to_decimal(x1_posix, profile_data$datetime)


  # Get corners of the parallelogram
  corners <- get_corners(x0_numeric, y0, x1_numeric, y1, slope)

  # Convert numeric x-values back to datetime for plotting
  corners_datetime <- lapply(corners, function(corner) {
    list(datetime = decimal_to_posixct(corner[1], profile_data$datetime),
         melatonin = corner[2])
  })


  # Create a dataframe for the parallelogram
  parallelogram_df <- do.call(rbind, lapply(corners_datetime, as.data.frame))

  # Add the parallelogram as a polygon to the plot
  plot <- plot +
    ggplot2::geom_polygon(
      data = parallelogram_df,
      ggplot2::aes(x = .data$datetime, y = .data$melatonin),
      fill = "red", alpha = 0.3
    )

  return(plot)
}

plot_profile <- function(profile_data, show_segments = TRUE, show_parallelogram = FALSE, pll_result = NULL) {
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

  # Add parallelogram overlay if show_parallelogram is TRUE
  if (show_parallelogram) {
    plot <- plot_parallelogram(plot, profile_data, pll_result)
  }

  # Return the plot object
  return(plot)
}

