# # last best version
# plot_base_segment <- function(profile_data) {
#   ggplot2::ggplot(profile_data, ggplot2::aes(x = .data$datetime, y = .data$melatonin)) +
#     # Plot a single continuous line for the full profile (mapped to "Full Profile")
#     ggplot2::geom_line(
#       ggplot2::aes(color = "Full Profile", group = 1),
#       size = 1.25
#     ) +
#     # Overlay points for the base segment (mapped to "Base Segment")
#     ggplot2::geom_line(
#       data = dplyr::filter(profile_data, .data$base == 1),
#       ggplot2::aes(x = .data$datetime, y = .data$melatonin, color = "Base Segment"),
#       size = 1.25
#     ) +
#     # Format the x-axis to show only time
#     ggplot2::scale_x_datetime(
#       labels = scales::date_format("%H:%M"),
#       date_breaks = "2 hours"
#     ) +
#     # Add plot labels
#     ggplot2::labs(
#       title = "Melatonin Profile with Base Segment Highlighted",
#       x = "Time",
#       y = "Melatonin Concentration"
#     ) +
#     ggplot2::theme_minimal() +
#     # Show the legend on the right
#     ggplot2::theme(legend.position = "right")
# }

# plot_roi <- function(plot, roi){
#  # add roi segment to plot
#   plot <- plot +
#       ggplot2::geom_line() +
#       ggplot2::geom_segment(
#         ggplot2::aes(x = roi$x_start, xend = roi$x_end,
#             y = -0.05, yend = -0.05),
#         color = "purple", size = 1.5
#       ) +
# # add roi rectangle to plot
#        ggplot2::geom_rect(
#          ggplot2::aes(xmin = roi$x_start, xmax = roi$x_end,
#              ymin = roi$y_min, ymax = roi$y_max),
#          fill = "purple", alpha = 0.01
#       )
#   return(plot)
# }

# add DLMO fit lines to plot
plot_fit<- function(plot, profile_data, dlmoFit){
  # convert between  posixct and decimal hours
  xstart_num = posixct_to_decimal(dplyr::filter(profile_data, base == 1)$datetime[1], profile_data$datetime)
  xend_num = posixct_to_decimal(tail(dplyr::filter(profile_data, ascending == 1)$datetime,n=1), profile_data$datetime)
  ipx_posix = decimal_to_posixct(dlmoFit$inflection_point$x, profile_data$datetime)

  # left fit line (by default always linear) #TODO insert warning if fit is not linear for left segment
  plot<- plot +
    ggplot2::geom_segment( # m * (x - poi_x) + poi_y
      x = dplyr::filter(profile_data, base == 1)$datetime[1],
      #TODO this works!!
      y = dlmoFit$left_params * (xstart_num - dlmoFit$inflection_point$x) + dlmoFit$inflection_point$y,
      #y = -.1119806 * (xstart_num - 20.58333) + 0.4,
      #TODO this works!!!
      xend = ipx_posix,
      #xend = decimal_to_posixct(20.58333, profile_data$datetime),
      # TODO this works!!
      yend = dlmoFit$inflection_point$y,
      #yend = 0.4,
      color = "pink",
      size = 1
    )

  # right fit line (either linear or parabolic)
  # check fit type (length 1 = linear, length 3 = parabolic)
  if(length(dlmoFit$right_params) == 1){
    # TODO
  }
  else{ # parabolic fit plot
    plot <- plot +
      ggplot2::geom_function(
        fun = function(x) {
          # Convert x (POSIXct) to decimal hours
          x_numeric <- posixct_to_decimal(x, profile_data$datetime)
          #x_numeric <- c(dlmoFit$inflection_point$x, xend_num)

          # Evaluate the polynomial function in decimal hours
          # TODO THIS WORKS!!!!
           y <- dlmoFit$right_params[1] * x_numeric^2 +
             dlmoFit$right_params[2] * x_numeric +
             dlmoFit$right_params[3]
          #y<- -3.378851 * x_numeric^2 + 159.363581 * x_numeric - 1848.303449
          return(y)
        },
        color = 'pink',
        size = 1,
        n = 1000,
        xlim = c(
          # TODO this works!
          decimal_to_posixct(dlmoFit$inflection_point$x, profile_data$datetime),
          #decimal_to_posixct(20.58333, profile_data$datetime),
          decimal_to_posixct(xend_num, profile_data$datetime)
        )
      )
  }
  ggplot2::scale_x_datetime(
    limits = c(decimal_to_posixct(dlmoFit$inflection_point$x, profile_data$datetime),
               decimal_to_posixct(xend_num, profile_data$datetime)),
    date_labels = "%H:%M", # Adjust date labels as needed
    date_breaks = "1 hour" # Adjust date breaks as needed
  )
return(plot)
}
# add DLMO inflection point to plot
plot_ip <- function(plot, profile_data, dlmoip){
  plot<- plot +
    ggplot2::geom_point(
      # TODO This works!!
      x = decimal_to_posixct(dlmoip$x, profile_data$datetime),
      #x = decimal_to_posixct(20.58333, profile_data$datetime),
      # TODO This works!!
      y = dlmoip$y,
      #y = 0.4,
      color = "darkorchid",
      size = 4,
      shape = 18
    )
  return(plot)
}

plot_roi <- function(plot, roi) {
  # Add ROI segment and rectangle to the plot
  plot <- plot +
    ggplot2::geom_segment(
      x = roi$x_start, xend = roi$x_end,
      y = -0.2, yend = -0.2,
      color = "#C77CFF", size = 1.5
    ) +
    ggplot2::geom_rect(
      xmin = roi$x_start, xmax = roi$x_end,
      ymin = roi$y_min, ymax = roi$y_max,
      fill = "#C77CFF", alpha = 0.01
    )
  return(plot)
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

plot_profile <- function(profile_data, show_segments = TRUE, show_parallelogram = FALSE, pll_result = NULL, show_roi = FALSE, roi = NULL, show_dlmoIP = TRUE, dlmoFit = NULL, show_fit = FALSE) {
  plot <- ggplot2::ggplot(profile_data, ggplot2::aes(x = .data$datetime, y = .data$melatonin)) +
    # Plot a single dotted line for the full profile (mapped to "Full Profile")
    ggplot2::geom_point(
      color = 'grey',
      size = 2
    ) + ggplot2::geom_line(color = 'grey', linetype = "dotted", size = 1)+
    # Format the x-axis to show only time
    ggplot2::scale_x_datetime(
      labels = scales::date_format("%H:%M"),
      date_breaks = "2 hours"
    ) +
    # Add plot labels
    ggplot2::labs(
      title = "Melatonin profile",
      x = "Local time [hh:mm]",
      y = "Melatonin concentration [pg/mL]"
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
      ggplot2::geom_point(
        data = dplyr::filter(profile_data, .data$base == 1),
        ggplot2::aes(x = .data$datetime, y = .data$melatonin, color = "Base Segment"),
        color = '#56B4E9',
        size = 2
      ) +
      # Overlay points for the ascending segment (mapped to "Ascending Segment")
      ggplot2::geom_point(
        data = dplyr::filter(profile_data, .data$ascending == 1),
        ggplot2::aes(x = .data$datetime, y = .data$melatonin, color = "Ascending Segment"),
        color = 'lightgreen',
        size = 2
      )
  }

  # Add parallelogram overlay if show_parallelogram is TRUE
  if (show_parallelogram) {
    plot <- plot_parallelogram(plot, profile_data, pll_result)
  }

  # Add region of interest overlay, if show_roi is TRUE
  if (show_roi){
    plot <- plot_roi(plot, roi)
  }

  # Add DLMO fit lines
  if (show_fit){
    plot<- plot_fit(plot, profile_data, dlmoFit)
  }

  # Add DLMO inflection point
  if (show_dlmoIP){
    plot<- plot_ip(plot,profile_data, dlmoFit$inflection_point)
  }
  # Return the plot object
  return(plot)
}

