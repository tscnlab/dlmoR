# # load data
# filename <- system.file("extdata/sample_melatonin_profile.csv", package = "dlmoR")
# sample_dlmo <- calculate_dlmo(file_path = filename, threshold = 2.3)
#
# # calculate time interval
# dh_time<-posixct_to_decimal(sample_dlmo$prof$datetime, sample_dlmo$prof$datetime)
# min_interval <- min(diff(sort(dh_time)))
#
# # set parameters
# time_sd <- min_interval/2 # SD of time noise in decimal hours, = 1/2 of smallest sampling interval
# mel_sd <- 0.3 # SD of melatonin noise in pg/mL for ELISA/RIA error
# n_iter <- 1000 # of noise realizations
#
# # store results
# dlmo_estimate_dh <- numeric(n_iter)
# dlmo_estimate_posix <- as.POSIXct(rep(NA, n_iter), origin = as.Date(sample_dlmo$prof$datetime[1]), tz = "UTC")
#
# # bootstrap loop
# set.seed(123)
# function to ensure
enforce_min_time_gap <- function(time_vec, min_gap = 1/60) {

  for (i in 2:length(time_vec)) {
    gap <- time_vec[i] - time_vec[i - 1]
    if (gap < min_gap) {
      time_vec[i] <- time_vec[i - 1] + min_gap
    }
  }
  return(time_vec)
}

#for (i in seq_len(n_iter)){
for (i in seq_len(n_iter)){
  # 1: add noise
  time_pert <- dh_time + rnorm(length(dh_time), mean = 0, sd = time_sd)
  mel_pert <- sample_dlmo$prof$melatonin + rnorm(length(sample_dlmo$prof$melatonin), mean = 0, sd = mel_sd)

  # 2: chronologically sort vectors
  ord <- order(time_pert)
  time_sorted <- time_pert[ord]
  mel_sorted <- mel_pert[ord]

  # 3: enforce time gaps between samples after perturbation
  time_sorted_gap<- enforce_min_time_gap(time_sorted, min_gap = 1 / 60)  # 1 minute in decimal hours

  # create df
  df <- tibble::tibble(
    datetime = decimal_to_posixct(time_sorted_gap, sample_dlmo$prof$datetime),
    melatonin = mel_sorted
  )

  #estimate dlmo
  iter_output <- calculate_dlmo(data = df)
  dlmo_estimate_posix[i] <- decimal_to_posixct(iter_output$ip$inflection_point_fine$x, sample_dlmo$prof$datetime)
  dlmo_estimate_dh[i] <- iter_output$ip$inflection_point_fine$x
}

summary(dlmo_estimate_dh)

# compute mean & confidence intervals

# decimal hours
dlmo_mean <- mean(dlmo_estimate_dh, na.rm = TRUE)
dlmo_ci <- quantile(dlmo_estimate_dh, probs = c(0.025, 0.975), na.rm = TRUE)

# posixct
dlmo_mean_posix <- decimal_to_posixct(dlmo_mean, sample_dlmo$prof$datetime)
dlmo_ci_posix <- decimal_to_posixct(dlmo_ci, sample_dlmo$prof$datetime)

# bootstrapped plots
# Prepare data and sample size

plot_df <- data.frame(dlmo = dlmo_estimate_posix)
actual_dlmo <- decimal_to_posixct(sample_dlmo$ip$inflection_point_fine$x, sample_dlmo$prof$datetime)


# Compute density and scale
density_vals <- density(as.numeric(plot_df$dlmo))
max_hist_count <- max(ggplot2::ggplot_build(
  ggplot2::ggplot(plot_df, ggplot2::aes(x = dlmo)) +
    ggplot2::geom_histogram(binwidth = 60 * 5))$data[[1]]$count)

y_text <- max_hist_count - 0.3  # position just above tallest bar

scaled_density <- data.frame(
  x = as.POSIXct(density_vals$x, origin = "1970-01-01", tz = "UTC"),
  y = density_vals$y * max_hist_count / max(density_vals$y)
)

# Prepare CI data for shading
dlmo_ci <- quantile(plot_df$dlmo, c(0.025, 0.975))
ci_df <- data.frame(
  xmin = dlmo_ci[1],
  xmax = dlmo_ci[2],
  ymin = 0,
  ymax = Inf,
  group = factor("95% CI", levels = c("DLMO estimate", "Bootstrap mean", "95% CI"))
)

# Prepare vertical lines
line_df <- data.frame(
  label = factor(c("DLMO estimate", "Bootstrap mean"),
                 levels = c("DLMO estimate", "Bootstrap mean", "95% CI")),
  x = c(actual_dlmo, mean(plot_df$dlmo))
)

# Final plot
ggplot2::ggplot(plot_df, ggplot2::aes(x = dlmo)) +
  # Histogram
  ggplot2::geom_histogram(
    binwidth = 60 * 5,
    fill = "lightblue",
    color = "black",
    boundary = 0
  ) +
  # Shaded 95% CI region
  ggplot2::geom_rect(
    data = ci_df,
    ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = group),
    alpha = 0.15,
    inherit.aes = FALSE
  ) +
  # Density curve
  ggplot2::geom_line(
    data = scaled_density,
    ggplot2::aes(x = x, y = y),
    color = "darkblue",
    linewidth = 1
  ) +
  # Vertical lines
  ggplot2::geom_vline(
    data = line_df,
    ggplot2::aes(xintercept = x, color = label, linetype = label),
    linewidth = 1
  ) +
  # Manual color/linetype legends
  ggplot2::scale_color_manual(
    values = c(
      "DLMO estimate" = "deeppink",
      "Bootstrap mean" = "blue"
    )
  ) +
  ggplot2::scale_linetype_manual(
    values = c(
      "DLMO estimate" = "dotdash",
      "Bootstrap mean" = "solid"
    )
  ) +
  ggplot2::scale_fill_manual(
    name = "Reference",
    values = c("95% CI" = "azure4")
  ) +
  ggplot2::scale_x_datetime(date_labels = "%H:%M", date_breaks = "1 hour") +
  ggplot2::scale_y_continuous(
    name = "Count",
    sec.axis = ggplot2::sec_axis(~ . * max(density_vals$y) / max_hist_count, name = "Density")
  ) +
  ggplot2::labs(
    title = "Bootstrapped DLMO Estimate: Histogram + Scaled Density",
    x = "DLMO Estimate [hh:mm]"
  ) +
  ggplot2::annotate(
    "text",
    x = actual_dlmo,
    y = y_text,
    label = paste0("DLMO: ", sample_dlmo$dlmo$fine$time),
    color = "deeppink",
    angle = 0,
    hjust = -5,
    vjust = 1,
    size = 3
  ) +
  ggplot2::annotate(
    "text",
    x = mean(plot_df$dlmo),
    y = y_text,
    label = paste0("Mean: ", format(dlmo_mean_posix, "%H:%M:%S")),
    color = "blue",
    angle = 0,
    hjust = -5,
    vjust = 2.5,
    size = 3
  )+
  ggplot2::theme_minimal() +
  ggplot2::guides(
    color = ggplot2::guide_legend(order = 1),
    linetype = ggplot2::guide_legend(order = 1),
    fill = ggplot2::guide_legend(order = 2)
  )+
  ggplot2::theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = ggplot2::element_blank(),
    legend.box = "horizontal",
    legend.margin = ggplot2::margin(t = 0)  # optional tweak to move closer to plot
  )


# ##
#
# df <- data.frame(dlmo = dlmo_estimate_posix)
# n <- nrow(df)
#
# # Plot
# ggplot2::ggplot(df, ggplot2::aes(x = dlmo)) +
#   ggplot2::geom_density(
#     ggplot2::aes(y = ggplot2::after_stat(density) * n),
#     fill = "lightblue", alpha = 0.6
#   ) +
#   ggplot2::geom_vline(xintercept = mean(df$dlmo), color = "blue", size = 1) +
#   ggplot2::geom_vline(
#     xintercept = quantile(df$dlmo, c(0.025, 0.975)),
#     color = "red", linetype = "dashed"
#   ) +
#   ggplot2::scale_x_datetime(date_labels = "%H:%M", date_breaks = "1 hour") +
#   ggplot2::labs(
#     title = "Bootstrapped DLMO Estimate Density",
#     x = "DLMO Estimate [hh:mm]",
#     y = "Estimated Count"
#   ) +
#   ggplot2::theme_minimal()
#
# ## histo
# df <- data.frame(dlmo = dlmo_estimate_posix)
#
# ggplot2::ggplot(df, ggplot2::aes(x = dlmo)) +
#   ggplot2::geom_histogram(
#     binwidth = 60 * 5,  # 5-minute bins (in seconds)
#     fill = "lightblue",
#     color = "black",
#     boundary = 0
#   ) +
#   ggplot2::geom_vline(
#     xintercept = mean(df$dlmo),
#     color = "blue",
#     size = 1
#   ) +
#   ggplot2::geom_vline(
#     xintercept = quantile(df$dlmo, c(0.025, 0.975)),
#     color = "red",
#     linetype = "dashed"
#   ) +
#   ggplot2::scale_x_datetime(
#     date_labels = "%H:%M",
#     date_breaks = "1 hour"
#   ) +
#   ggplot2::labs(
#     title = "Bootstrapped DLMO Estimate Histogram",
#     x = "DLMO Estimate [hh:mm]",
#     y = "Count"
#   ) +
#   ggplot2::theme_minimal()
#
# ## histo + density
# df <- data.frame(dlmo = dlmo_estimate_posix)
# n <- nrow(df)
#
# ggplot2::ggplot(df, ggplot2::aes(x = dlmo)) +
#   # Histogram in counts
#   ggplot2::geom_histogram(
#     binwidth = 60 * 5,  # 5-minute bins
#     fill = "lightblue",
#     color = "black",
#     boundary = 0
#   ) +
#   # Density overlaid and scaled to match count scale
#   ggplot2::geom_density(
#     ggplot2::aes(y = ..density.. * n),
#     color = "darkblue",
#     linewidth = 1
#   ) +
#   # Mean DLMO line
#   ggplot2::geom_vline(
#     xintercept = mean(df$dlmo),
#     color = "blue",
#     size = 1
#   ) +
#   # 95% CI lines
#   ggplot2::geom_vline(
#     xintercept = quantile(df$dlmo, c(0.025, 0.975)),
#     color = "red",
#     linetype = "dashed"
#   ) +
#   # Time axis formatting
#   ggplot2::scale_x_datetime(
#     date_labels = "%H:%M",
#     date_breaks = "1 hour"
#   ) +
#   # Labels
#   ggplot2::labs(
#     title = "Bootstrapped DLMO Estimate Histogram + Density",
#     x = "DLMO Estimate [hh:mm]",
#     y = "Count"
#   ) +
#   ggplot2::theme_minimal()
#
# ## try 2
#
# df <- data.frame(dlmo = dlmo_estimate_posix)
# actual_dlmo <- decimal_to_posixct(sample_dlmo$ip$inflection_point_fine$x, sample_dlmo$prof$datetime)
# # Manually compute density
# density_vals <- density(as.numeric(df$dlmo))  # POSIXct auto-converts to seconds
# max_hist_count <- max(ggplot2::ggplot_build(
#   ggplot2::ggplot(df, ggplot2::aes(x = dlmo)) +
#     ggplot2::geom_histogram(binwidth = 60 * 5))$data[[1]]$count)
#
# scaled_density <- data.frame(
#   x = as.POSIXct(density_vals$x, origin = "1970-01-01", tz = "UTC"),
#   y = density_vals$y * max_hist_count / max(density_vals$y)  # Scale to histogram peak
# )
#
# # Plot histogram + scaled density
# ggplot2::ggplot(df, ggplot2::aes(x = dlmo)) +
#   ggplot2::geom_histogram(
#     binwidth = 60 * 5,
#     fill = "lightblue",
#     color = "black",
#     boundary = 0
#   ) +
#   ggplot2::geom_line(
#     data = scaled_density,
#     ggplot2::aes(x = x, y = y),
#     color = "darkblue",
#     linewidth = 1
#   ) +
#   ggplot2::geom_vline(xintercept = mean(df$dlmo), color = "blue", size = 1) +
#   ggplot2::geom_vline(xintercept = quantile(df$dlmo, c(0.025, 0.975)), color = "red", linetype = "dashed") +
#
#   # Add this to your ggplot:
#   ggplot2::geom_vline(
#     xintercept = actual_dlmo,
#     color = "magenta4",
#     linetype = "dotdash",
#     linewidth = 1
#   )+
#
#   ggplot2::scale_x_datetime(date_labels = "%H:%M", date_breaks = "1 hour") +
#   ggplot2::labs(
#     title = "Bootstrapped DLMO Estimate: Histogram + Scaled Density",
#     x = "DLMO Estimate [hh:mm]",
#     y = "Count"
#   ) +
#   ggplot2::annotate(
#     "text",
#     x = actual_dlmo,
#     y = Inf,
#     label = "Original DLMO",
#     vjust = 2,
#     angle = 90,
#     color = "magenta4",
#     size = 3
#   )+
#   ggplot2::theme_minimal()
#
# ##
# df <- data.frame(dlmo = dlmo_estimate_posix)
#
# density_vals <- density(as.numeric(df$dlmo))  # POSIXct → numeric (seconds)
# max_hist_count <- max(ggplot2::ggplot_build(
#   ggplot2::ggplot(df, ggplot2::aes(x = dlmo)) +
#     ggplot2::geom_histogram(binwidth = 60 * 5))$data[[1]]$count)
#
# scaled_density <- data.frame(
#   x = as.POSIXct(density_vals$x, origin = sample_dlmo$prof$datetime, tz = "UTC"),
#   y = density_vals$y * max_hist_count / max(density_vals$y)
# )
#
# ggplot2::ggplot(df, ggplot2::aes(x = dlmo)) +
#   ggplot2::geom_histogram(
#     binwidth = 60 * 5,
#     fill = "lightblue",
#     color = "black",
#     boundary = 0
#   ) +
#   ggplot2::geom_line(
#     data = scaled_density,
#     ggplot2::aes(x = x, y = y),
#     color = "darkblue",
#     linewidth = 1
#   ) +
#   ggplot2::geom_vline(xintercept = mean(df$dlmo), color = "blue", size = 1) +
#   ggplot2::geom_vline(xintercept = quantile(df$dlmo, c(0.025, 0.975)), color = "red", linetype = "dashed") +
#   ggplot2::scale_x_datetime(date_labels = "%H:%M", date_breaks = "1 hour") +
#   ggplot2::scale_y_continuous(
#     name = "Count",
#     sec.axis = ggplot2::sec_axis(~ . * max(density_vals$y) / max_hist_count, name = "Density")
#   ) +
#   ggplot2::labs(
#     title = "Bootstrapped DLMO Estimate: Histogram + Scaled Density",
#     x = "DLMO Estimate [hh:mm]"
#   ) +
#   ggplot2::theme_minimal()
#
# ##
# df <- data.frame(dlmo = dlmo_estimate_posix)
# actual_dlmo <- decimal_to_posixct(sample_dlmo$ip$inflection_point_fine$x, sample_dlmo$prof$datetime)
#
# # Manually compute density
# density_vals <- density(as.numeric(df$dlmo))
# max_hist_count <- max(ggplot2::ggplot_build(
#   ggplot2::ggplot(df, ggplot2::aes(x = dlmo)) +
#     ggplot2::geom_histogram(binwidth = 60 * 5))$data[[1]]$count)
#
# scaled_density <- data.frame(
#   x = as.POSIXct(density_vals$x, origin = sample_dlmo$prof$datetime, tz = "UTC"),
#   y = density_vals$y * max_hist_count / max(density_vals$y)
# )
#
# # Key summary values
# dlmo_mean <- mean(df$dlmo)
# dlmo_ci <- quantile(df$dlmo, c(0.025, 0.975))
#
# # Final plot
# ggplot2::ggplot(df, ggplot2::aes(x = dlmo)) +
#   ggplot2::geom_histogram(
#     binwidth = 60 * 5,
#     fill = "lightblue",
#     color = "black",
#     boundary = 0
#   ) +
#   ggplot2::geom_line(
#     data = scaled_density,
#     ggplot2::aes(x = x, y = y),
#     color = "darkblue",
#     linewidth = 1
#   ) +
#   # DLMO estimate (original)
#   ggplot2::geom_vline(xintercept = actual_dlmo, color = "magenta4", linetype = "dotdash", linewidth = 1) +
#   ggplot2::annotate(
#     "text",
#     x = actual_dlmo,
#     y = Inf,
#     label = "DLMO estimate",
#     vjust = 2, angle = 90,
#     color = "magenta4", size = 3
#   ) +
#   # Mean of bootstrapped estimates
#   ggplot2::geom_vline(xintercept = dlmo_mean, color = "blue", size = 1) +
#   ggplot2::annotate(
#     "text",
#     x = dlmo_mean,
#     y = Inf,
#     label = "Bootstrap mean",
#     vjust = 2, angle = 90,
#     color = "blue", size = 3
#   ) +
#   # Confidence interval lines
#   ggplot2::geom_vline(xintercept = dlmo_ci, color = "red", linetype = "dashed") +
#   ggplot2::annotate(
#     "text",
#     x = dlmo_ci[1],
#     y = Inf,
#     label = "2.5%",
#     vjust = 2, angle = 90,
#     color = "red", size = 3
#   ) +
#   ggplot2::annotate(
#     "text",
#     x = dlmo_ci[2],
#     y = Inf,
#     label = "97.5%",
#     vjust = 2, angle = 90,
#     color = "red", size = 3
#   ) +
#   ggplot2::scale_x_datetime(date_labels = "%H:%M", date_breaks = "1 hour") +
#   ggplot2::scale_y_continuous(
#     name = "Count",
#     sec.axis = ggplot2::sec_axis(~ . * max(density_vals$y) / max_hist_count, name = "Density")
#   ) +
#   ggplot2::labs(
#     title = "Bootstrapped DLMO Estimate: Histogram + Scaled Density",
#     x = "DLMO Estimate [hh:mm]"
#   ) +
#   ggplot2::theme_minimal()
#
# ##
# df <- data.frame(dlmo = dlmo_estimate_posix)
# actual_dlmo <- decimal_to_posixct(sample_dlmo$ip$inflection_point_fine$x, sample_dlmo$prof$datetime)
#
# # Density calculation (stay numeric here)
# density_vals <- density(as.numeric(df$dlmo))
# max_hist_count <- max(ggplot2::ggplot_build(
#   ggplot2::ggplot(df, ggplot2::aes(x = dlmo)) +
#     ggplot2::geom_histogram(binwidth = 60 * 5))$data[[1]]$count)
#
# scaled_density <- data.frame(
#   x = as.POSIXct(density_vals$x, origin = "1970-01-01", tz = "UTC"),
#   y = density_vals$y * max_hist_count / max(density_vals$y)
# )
#
# # Final plot
# ggplot2::ggplot(df, ggplot2::aes(x = dlmo)) +
#   ggplot2::geom_histogram(
#     binwidth = 60 * 5,
#     fill = "lightblue",
#     color = "black",
#     boundary = 0
#   ) +
#   ggplot2::geom_line(
#     data = scaled_density,
#     ggplot2::aes(x = x, y = y),
#     color = "darkblue",
#     linewidth = 1
#   ) +
#   ggplot2::geom_vline(xintercept = mean(df$dlmo), color = "blue", size = 1) +
#   ggplot2::geom_vline(xintercept = quantile(df$dlmo, c(0.025, 0.975)),
#                       color = "red", linetype = "dashed") +
#   ggplot2::geom_vline(xintercept = actual_dlmo,
#                       color = "magenta4", linetype = "dotdash", linewidth = 1) +
#
#   # Annotations
#   ggplot2::annotate("text", x = actual_dlmo, y = Inf, label = "DLMO estimate",
#                     vjust = -0.5, angle = 90, color = "magenta4", size = 3) +
#   ggplot2::annotate("text", x = mean(df$dlmo), y = Inf, label = "Mean",
#                     vjust = -0.5, angle = 90, color = "blue", size = 3) +
#   ggplot2::annotate("text", x = quantile(df$dlmo, 0.025), y = Inf, label = "2.5%",
#                     vjust = -0.5, angle = 90, color = "red", size = 3) +
#   ggplot2::annotate("text", x = quantile(df$dlmo, 0.975), y = Inf, label = "97.5%",
#                     vjust = -0.5, angle = 90, color = "red", size = 3) +
#
#   # Axes and theme
#   ggplot2::scale_x_datetime(date_labels = "%H:%M", date_breaks = "1 hour") +
#   ggplot2::scale_y_continuous(
#     name = "Count",
#     sec.axis = ggplot2::sec_axis(~ . * max(density_vals$y) / max_hist_count, name = "Density")
#   ) +
#   ggplot2::labs(
#     title = "Bootstrapped DLMO Estimate: Histogram + Scaled Density",
#     x = "DLMO Estimate [hh:mm]"
#   ) +
#   ggplot2::theme_minimal()
#
# ##
# df <- data.frame(dlmo = dlmo_estimate_posix)
# actual_dlmo <- decimal_to_posixct(sample_dlmo$ip$inflection_point_fine$x, sample_dlmo$prof$datetime)
#
# # Manually compute density
# density_vals <- density(as.numeric(df$dlmo))
# max_hist_count <- max(ggplot2::ggplot_build(
#   ggplot2::ggplot(df, ggplot2::aes(x = dlmo)) +
#     ggplot2::geom_histogram(binwidth = 60 * 5))$data[[1]]$count)
#
# scaled_density <- data.frame(
#   x = as.POSIXct(density_vals$x, origin = "1970-01-01", tz = "UTC"),
#   y = density_vals$y * max_hist_count / max(density_vals$y)
# )
#
# # Prepare summary values for line plotting
# line_df <- data.frame(
#   label = factor(c("DLMO estimate", "Bootstrap mean", "CI lower", "CI upper"),
#                  levels = c("DLMO estimate", "Bootstrap mean", "CI lower", "CI upper")),
#   x = c(actual_dlmo,
#         mean(df$dlmo),
#         quantile(df$dlmo, 0.025),
#         quantile(df$dlmo, 0.975))
# )
#
# # Plot with legend
# ggplot2::ggplot(df, ggplot2::aes(x = dlmo)) +
#   ggplot2::geom_histogram(
#     binwidth = 60 * 5,
#     fill = "lightblue",
#     color = "black",
#     boundary = 0
#   ) +
#   ggplot2::geom_line(
#     data = scaled_density,
#     ggplot2::aes(x = x, y = y),
#     color = "darkblue",
#     linewidth = 1
#   ) +
#   ggplot2::geom_vline(
#     data = line_df,
#     ggplot2::aes(xintercept = x, color = label, linetype = label),
#     linewidth = 1
#   ) +
#   ggplot2::scale_color_manual(
#     values = c("DLMO estimate" = "deeppink4",
#                "Bootstrap mean" = "blue",
#                "CI lower" = "red",
#                "CI upper" = "red")
#   ) +
#   ggplot2::scale_linetype_manual(
#     values = c("DLMO estimate" = "dotdash",
#                "Bootstrap mean" = "solid",
#                "CI lower" = "dashed",
#                "CI upper" = "dashed")
#   ) +
#   ggplot2::scale_x_datetime(date_labels = "%H:%M", date_breaks = "1 hour") +
#   ggplot2::scale_y_continuous(
#     name = "Count",
#     sec.axis = ggplot2::sec_axis(~ . * max(density_vals$y) / max_hist_count, name = "Density")
#   ) +
#   ggplot2::labs(
#     title = "Bootstrapped DLMO Estimate: Histogram + Scaled Density",
#     x = "DLMO Estimate [hh:mm]",
#     color = "Reference",  # Legend title
#     linetype = "Reference"
#   ) +
#   ggplot2::theme_minimal() +
#   ggplot2::theme(legend.position = "right")
