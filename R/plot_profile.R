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
#     # Show the legend on the ascending
#     ggplot2::theme(legend.position = "ascending")
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

  # base fit line (by default always linear) #TODO insert warning if fit is not linear for base segment
  plot<- plot +
    ggplot2::geom_segment( # m * (x - poi_x) + poi_y
      x = dplyr::filter(profile_data, base == 1)$datetime[1],
      #TODO this works!!
      y = dlmoFit$base_params * (xstart_num - dlmoFit$inflection_point$x) + dlmoFit$inflection_point$y,
      #y = -.1119806 * (xstart_num - 20.58333) + 0.4,
      #TODO this works!!!
      xend = ipx_posix,
      #xend = decimal_to_posixct(20.58333, profile_data$datetime),
      # TODO this works!!
      yend = dlmoFit$inflection_point$y,
      #yend = 0.4,
      color = "deeppink",
      size = 1
    )

  # ascending fit line (either linear or parabolic)
  # check fit type (length 1 = linear, length 3 = parabolic)
  if(length(dlmoFit$ascending_params) == 1){
    plot<- plot +
      ggplot2::geom_segment( # m * (poi_x - x) + poi_y
        #ascending_data = dplyr::filter(profile_data, ascending == 1)$datetime,
        x = ipx_posix,
        #TODO this works!!
        y = dlmoFit$inflection_point$y,
        #y = dlmoFit$ascending_params * (dlmoFit$inflection_point$x - xend_num) + dlmoFit$inflection_point$y,
        #y = -.1119806 * (xstart_num - 20.58333) + 0.4,
        #TODO this works!!!
        xend = dplyr::filter(profile_data, ascending == 1)$datetime[length(dplyr::filter(profile_data, ascending == 1)$datetime)],
        #xend = decimal_to_posixct(20.58333, profile_data$datetime),
        # TODO this works!!
        #yend = dlmoFit$inflection_point$y,
        yend = dplyr::filter(profile_data, ascending == 1)$melatonin[length(dplyr::filter(profile_data, ascending == 1)$melatonin)],
        #yend = 0.4,
        color = "deeppink",
        size = 1
      )
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
           y <- dlmoFit$ascending_params$a * x_numeric^2 +
             dlmoFit$ascending_params$b * x_numeric +
             dlmoFit$ascending_params$c
          #y<- -3.378851 * x_numeric^2 + 159.363581 * x_numeric - 1848.303449
          return(y)
        },
        color = 'deeppink',
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
      color = "deeppink",
      fill = "deeppink",
      size = 3,
      shape = 23
    )
  return(plot)
}

plot_roi <- function(plot, roi, roi_line_only) {
  # Add ROI segment and rectangle to the plot
  plot <- plot +
    ggplot2::geom_segment(
      x = roi$x_start, xend = roi$x_end,
      y = -0.2, yend = -0.2,
      color = "orchid3", size = 1.5
    ) +
    if(!roi_line_only){
    ggplot2::geom_rect(
      xmin = roi$x_start, xmax = roi$x_end,
      ymin = roi$y_min, ymax = roi$y_max,
      fill = "orchid3", alpha = 0.01
    )}
  return(plot)
}

plot_roi_heatmap <- function(plot, data = NULL, roi_grid_big = NULL, roi_grid_small = NULL, residuals_big = NULL, residuals_small = NULL, show_roi_big = FALSE, show_roi_small = FALSE){ #TODO make this cleaner
  if(show_roi_big & !is.null(roi_grid_big)){
print("allo!")
  # convert grid x coordinate to datetime
    dt_roi_grid_big <- data.frame(
      x = decimal_to_posixct(roi_grid_big$x, data$datetime),
      y = roi_grid_big$y
    )
  #print(dt_roi_grid_big)
  plot <- plot + ggplot2::geom_point(data=dt_roi_grid_big, ggplot2::aes(x = x, y = y, color = log(residuals_big))) +
    ggplot2::scale_color_gradient(low = "deeppink", high = "cyan")  # Color gradient from blue to red
    #ggplot2::scale_color_gradient(low = "orchid1", high = "orchid4")  # Color gradient from blue to red
  } else if(show_roi_small & !is.null(roi_grid_small)){
print("hi!")
    # convert grid x coordinate to datetime
    dt_roi_grid_small <- data.frame(
      x = decimal_to_posixct(roi_grid_small$x, data$datetime),
      y = roi_grid_small$y
    )

    plot <- plot + ggplot2::geom_point(data=dt_roi_grid_small, ggplot2::aes(x = x, y = y, color = log(residuals_small)), size = 0.5) +
      ggplot2::scale_color_gradient(low = "deeppink",
                                    high = "cyan",
                                    guide = "colorbar")  # Color gradient from blue to red
  }
  }


plot_parallelogram <- function(plot, profile_data, pll_result) {
  if (is.null(pll_result)) {
    stop("pll_result must be provided to plot the parallelogram.")
  }

  # Extract optimized parameters
  x0_posix <- pll_result$pll_datetime_0
  x1_posix <- pll_result$pll_datetime_1
  slope <- pll_result$pll_slope

  # # Filter for ascending data
  # profile_data_ascending <- profile_data %>% dplyr::filter(.data$ascending == 1)
  # y0 <- min(profile_data_ascending$melatonin)
  # y1 <- max(profile_data_ascending$melatonin)
  #
  # # Convert datetime to numeric for parallelogram calculations
  # x0_numeric <- posixct_to_decimal(x0_posix, profile_data$datetime)
  # x1_numeric <- posixct_to_decimal(x1_posix, profile_data$datetime)
  #
  #
  # # Get corners of the parallelogram
  # corners <- get_corners(x0_numeric, y0, x1_numeric, y1, slope)
    corners <- pll_result$corners
  # Convert numeric x-values back to datetime for plotting
  corners_datetime <- lapply(corners, function(corner) {
    list(datetime = decimal_to_posixct(corner[1], profile_data$datetime),
         melatonin = corner[2])
  })


  # Create a dataframe for the parallelogram
  parallelogram_df <- do.call(rbind, lapply(corners_datetime, as.data.frame))

  # TODO commented out 13.12.2024
  # Add the parallelogram as a polygon to the plot
  # plot <- plot +
  #   ggplot2::geom_polygon(
  #     data = parallelogram_df,
  #     ggplot2::aes(x = .data$datetime, y = .data$melatonin),
  #     fill = "red", alpha = 0.3
  #   )
  # Extract diagonal points
  diagonal_1 <- parallelogram_df[c(1, 3), ]  # Connect corner 1 and 3
  diagonal_2 <- parallelogram_df[c(2, 4), ]  # Connect corner 2 and 4

  # Add the parallelogram as a polygon to the plot
  plot <- plot +
    ggplot2::geom_polygon(
      data = parallelogram_df,
      ggplot2::aes(x = .data$datetime, y = .data$melatonin),
      fill = "red", alpha = 0.1
    ) +
    # Add diagonals as lines
    ggplot2::geom_line(
      data = diagonal_1,
      ggplot2::aes(x = .data$datetime, y = .data$melatonin),
      color = "red", linetype = "dashed", size = 0.5
    ) +
    ggplot2::geom_line(
      data = diagonal_2,
      ggplot2::aes(x = .data$datetime, y = .data$melatonin),
      color = "red", linetype = "dashed", size = 0.5
    )


  return(plot)
}

plot_profile <- function(profile_data, show_threshold = TRUE, threshold = 2.3, show_segments = TRUE, show_parallelogram = FALSE, pll_result = NULL, show_roi = FALSE, roi_line_only = TRUE, roi = NULL, show_dlmoIP = TRUE, dlmoFit = NULL, show_fit = FALSE, show_roi_heatmap = FALSE, show_roi_small = FALSE, show_roi_big = FALSE) {
  # Define the title conditionally
  plot_title <- if (!is.null(dlmoFit)) {
    dlmo_time <- hms::as_hms(decimal_to_posixct(dlmoFit$inflection_point$x, profile_data$datetime))

    paste("Melatonin profile\nDLMO time:", as.character(dlmo_time))
  } else {
    "Melatonin profile"
  }

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
      title = plot_title,
      x = "Local time [hh:mm]",
      y = "Melatonin concentration [pg/mL]"
    ) +
    ggplot2::theme_minimal() +
    # Show the legend on the ascending
    ggplot2::theme(legend.position = "ascending") +
    # Customize line types in the legend
    ggplot2::scale_linetype_manual(values = c("Full Profile" = "dotted"))




  # Add parallelogram overlay if show_parallelogram is TRUE
  if (show_parallelogram) {
    plot <- plot_parallelogram(plot, profile_data, pll_result)
  }

  # Add region of interest overlay, if show_roi is TRUE
  if (show_roi){
    plot <- plot_roi(plot, roi, roi_line_only)
  }

  # Add ROI heatmap, if show_roi_heatmap is TRUE

  if (show_roi_heatmap){
    plot <- plot_roi_heatmap(plot, data = profile_data, roi_grid_big = dlmoFit$grid_big, roi_grid_small = dlmoFit$grid_small, residuals_big = dlmoFit$res_big,residuals_small = dlmoFit$res_small, show_roi_small = show_roi_small, show_roi_big = show_roi_big)
  }

  # Plot a single dotted line for the full profile (mapped to "Full Profile")
  plot<- plot + ggplot2::geom_line(color = 'grey', linetype = "dotted", size = 1)

  # Add DLMO fit lines
  if (show_fit){
    plot<- plot_fit(plot, profile_data, dlmoFit)
  }

  # Add DLMO inflection point
  if (show_dlmoIP){
    plot<- plot_ip(plot,profile_data, dlmoFit$inflection_point)
  }

  # Add threshold line
  if (show_threshold){
    plot<- plot + ggplot2::geom_hline(ggplot2::aes(yintercept = threshold), color = "burlywood3", size = 1)
  }

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
      )+
      if("intermediate"%in%colnames(profile_data)){
        # Overlay points for the intermediate segment (mapped to "Intermediate Segment")
        ggplot2::geom_point(
          data = dplyr::filter(profile_data, .data$intermediate == 1),
          ggplot2::aes(x = .data$datetime, y = .data$melatonin, color = "Intermediate Segment"),
          color = 'darkgoldenrod1',
          size = 2
        )
      }
  }




  # Return the plot object
  return(plot)
}

