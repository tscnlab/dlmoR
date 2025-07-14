# Load required libraries
library(tibble)
library(ggplot2)

# Calculate DLMO and extract time
filename <- system.file("extdata/sample_melatonin_profile.csv", package = "dlmoR")
sample_dlmo <- calculate_dlmo(file_path = filename, threshold = 2.3)

# Convert time to decimal hours
dh_time <- posixct_to_decimal(sample_dlmo$prof$datetime, sample_dlmo$prof$datetime)
min_interval <- min(diff(sort(dh_time)))

# Set parameters
time_sd <- min_interval / 2  # time noise (decimal hours)
mel_sd <- 0.3                # melatonin noise (pg/mL)
n_iter <- 100               # bootstrap iterations

# Output storage
dlmo_estimate_dh <- numeric(n_iter)
dlmo_estimate_posix <- as.POSIXct(rep(NA, n_iter),
                                  origin = as.Date(sample_dlmo$prof$datetime[1]), tz = "UTC")

# Define minimum time gap function
enforce_min_time_gap <- function(time_vec, min_gap = 1 / 60) {
  for (i in 2:length(time_vec)) {
    gap <- time_vec[i] - time_vec[i - 1]
    if (gap < min_gap) {
      time_vec[i] <- time_vec[i - 1] + min_gap
    }
  }
  return(time_vec)
}

# PCHIP-style interpolator
interp_fun <- splinefun(
  x = dh_time,
  y = sample_dlmo$prof$melatonin,
  method = "monoH.FC"
)

# Bootstrap loop
set.seed(123)
for (i in seq_len(n_iter)) {
  time_pert <- dh_time + rnorm(length(dh_time), mean = 0, sd = time_sd)
  time_sorted <- enforce_min_time_gap(sort(time_pert), min_gap = 1 / 60)

  mel_interp <- interp_fun(time_sorted)
  mel_pert <- mel_interp + rnorm(length(mel_interp), mean = 0, sd = mel_sd)

  df <- tibble(
    datetime = decimal_to_posixct(time_sorted, sample_dlmo$prof$datetime),
    melatonin = mel_pert
  )

  iter_output <- calculate_dlmo(data = df)
  dlmo_estimate_posix[i] <- decimal_to_posixct(iter_output$ip$inflection_point_fine$x, sample_dlmo$prof$datetime)
  dlmo_estimate_dh[i] <- iter_output$ip$inflection_point_fine$x
}

# Summary statistics
dlmo_mean <- mean(dlmo_estimate_dh, na.rm = TRUE)
dlmo_ci <- quantile(dlmo_estimate_dh, probs = c(0.025, 0.975), na.rm = TRUE)
dlmo_mean_posix <- decimal_to_posixct(dlmo_mean, sample_dlmo$prof$datetime)
dlmo_ci_posix <- decimal_to_posixct(dlmo_ci, sample_dlmo$prof$datetime)

# --- Plot: Visual validation of interpolation ---
# Create time grid
time_grid <- seq(min(dh_time), max(dh_time), length.out = 500)
mel_interp_grid <- interp_fun(time_grid)

# Create one bootstrap sample
set.seed(42)
time_pert <- dh_time + rnorm(length(dh_time), mean = 0, sd = time_sd)
time_sorted <- enforce_min_time_gap(sort(time_pert), min_gap = 1 / 60)
mel_interp <- interp_fun(time_sorted)
mel_pert <- mel_interp + rnorm(length(mel_interp), mean = 0, sd = mel_sd)

# Convert to datetime
datetime_grid <- decimal_to_posixct(time_grid, sample_dlmo$prof$datetime)
datetime_sorted <- decimal_to_posixct(time_sorted, sample_dlmo$prof$datetime)

# Validation plot
# Prepare data for legend-enabled plotting

# Original profile
original_df <- sample_dlmo$prof %>%
  dplyr::mutate(source = "Original measurements")

# Interpolated curve from dense time grid
interp_df <- tibble::tibble(
  datetime = decimal_to_posixct(time_grid, sample_dlmo$prof$datetime),
  melatonin = mel_interp_grid,
  source = "Interpolated curve"
)

# One bootstrap perturbation
perturbed_df <- tibble::tibble(
  datetime = decimal_to_posixct(time_sorted, sample_dlmo$prof$datetime),
  melatonin = mel_pert,
  source = "Perturbed + noisy sample"
)

# Combine all
plot_df <- dplyr::bind_rows(original_df, interp_df, perturbed_df)

# Create the plot
ggplot(plot_df, aes(x = datetime, y = melatonin, color = source)) +
  # Interpolated curve (line only)
  geom_line(
    data = dplyr::filter(plot_df, source == "Interpolated curve"),
    linewidth = 1
  ) +
  # Points: original + perturbed
  geom_point(
    data = dplyr::filter(plot_df, source != "Interpolated curve"),
    size = 2, alpha = 0.8
  ) +
  scale_color_manual(
    name = "Data Source",
    values = c(
      "Original measurements" = "black",
      "Interpolated curve" = "blue",
      "Perturbed + noisy sample" = "red"
    )
  ) +
  labs(
    title = "Validation of PCHIP-like Interpolation for Bootstrap Sampling",
    x = "Time",
    y = "Melatonin (pg/mL)"
  ) +
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "1 hour") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

# Only proceed if you have enough bootstrap samples
if (n_iter > 1) {

  # Histogram + density plot of bootstrapped DLMO estimates
  plot_df <- data.frame(dlmo = dlmo_estimate_posix)
  actual_dlmo <- decimal_to_posixct(sample_dlmo$ip$inflection_point_fine$x, sample_dlmo$prof$datetime)

  # Density and scaling
  density_vals <- density(as.numeric(plot_df$dlmo), na.rm = TRUE)
  max_hist_count <- max(ggplot_build(
    ggplot(plot_df, aes(x = dlmo)) +
      geom_histogram(binwidth = 60 * 5))$data[[1]]$count)

  y_text <- max_hist_count - 0.3  # For annotation

  scaled_density <- data.frame(
    x = as.POSIXct(density_vals$x, origin = "1970-01-01", tz = "UTC"),
    y = density_vals$y * max_hist_count / max(density_vals$y)
  )

  # CI shading region
  ci_df <- data.frame(
    xmin = dlmo_ci_posix[1],
    xmax = dlmo_ci_posix[2],
    ymin = 0,
    ymax = Inf,
    group = factor("95% CI", levels = c("DLMO estimate", "Bootstrap mean", "95% CI"))
  )

  # Vertical lines
  line_df <- data.frame(
    label = factor(c("DLMO estimate", "Bootstrap mean"),
                   levels = c("DLMO estimate", "Bootstrap mean", "95% CI")),
    x = c(actual_dlmo, dlmo_mean_posix)
  )

  # Final plot
  ggplot(plot_df, aes(x = dlmo)) +
    geom_histogram(
      binwidth = 60 * 5,
      fill = "lightblue",
      color = "black",
      boundary = 0
    ) +
    geom_rect(
      data = ci_df,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = group),
      alpha = 0.15,
      inherit.aes = FALSE
    ) +
    geom_line(
      data = scaled_density,
      aes(x = x, y = y),
      color = "darkblue",
      linewidth = 1
    ) +
    geom_vline(
      data = line_df,
      aes(xintercept = x, color = label, linetype = label),
      linewidth = 1
    ) +
    scale_color_manual(
      values = c("DLMO estimate" = "deeppink", "Bootstrap mean" = "blue")
    ) +
    scale_linetype_manual(
      values = c("DLMO estimate" = "dotdash", "Bootstrap mean" = "solid")
    ) +
    scale_fill_manual(
      name = "Reference",
      values = c("95% CI" = "azure4")
    ) +
    scale_x_datetime(date_labels = "%H:%M", date_breaks = "1 hour") +
    scale_y_continuous(
      name = "Count",
      sec.axis = sec_axis(~ . * max(density_vals$y) / max_hist_count, name = "Density")
    ) +
    labs(
      title = "Bootstrapped DLMO Estimate: Histogram + Scaled Density",
      x = "DLMO Estimate [hh:mm]"
    ) +
    annotate(
      "text",
      x = actual_dlmo,
      y = y_text,
      label = paste0("DLMO: ", sample_dlmo$dlmo$fine$time),
      color = "deeppink",
      angle = 0,
      hjust = -0.1,
      vjust = 1,
      size = 3
    ) +
    annotate(
      "text",
      x = dlmo_mean_posix,
      y = y_text,
      label = paste0("Mean: ", format(dlmo_mean_posix, "%H:%M:%S")),
      color = "blue",
      angle = 0,
      hjust = -0.1,
      vjust = 2.5,
      size = 3
    ) +
    theme_minimal() +
    guides(
      color = guide_legend(order = 1),
      linetype = guide_legend(order = 1),
      fill = guide_legend(order = 2)
    ) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_blank(),
      legend.box = "horizontal",
      legend.margin = margin(t = 0)
    )
}

