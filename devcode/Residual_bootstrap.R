# --- Prerequisites ---
library(ggplot2)
library(dplyr)
library(tibble)

# Load sample DLMO profile
filename <- system.file("extdata/sample_melatonin_profile.csv", package = "dlmoR")
sample_dlmo <- calculate_dlmo(file_path = filename, threshold = 2.3)

profile_data <- sample_dlmo$prof
dlmoFit <- sample_dlmo$ip

# Convert datetimes to decimal hours
decimal_times <- posixct_to_decimal(profile_data$datetime, profile_data$datetime)
ip_x <- dlmoFit$inflection_point_fine$x
ip_y <- dlmoFit$inflection_point_fine$y

# Split data based on inflection point
# base_data <- profile_data[decimal_times <= ip_x, ]
# asc_data  <- profile_data[decimal_times >= ip_x, ]
#
# base_x <- posixct_to_decimal(base_data$datetime, profile_data$datetime)
# asc_x  <- posixct_to_decimal(asc_data$datetime, profile_data$datetime)


# Find start (first base point) and end (last ascending point)
first_base_idx <- which(profile_data$base == 1)[1]
last_asc_idx   <- tail(which(profile_data$ascending == 1), n = 1)

# Subset the data between those two points (inclusive)
middle_range <- profile_data[first_base_idx:last_asc_idx, ]
middle_times <- decimal_times[first_base_idx:last_asc_idx]

# Split based on inflection point
base_data <- middle_range[middle_times <= ip_x, ]
asc_data  <- middle_range[middle_times >= ip_x, ]

# base_x <- middle_times[middle_times <= ip_x]
# asc_x  <- middle_times[middle_times >= ip_x]
base_x <- posixct_to_decimal(base_data$datetime, profile_data$datetime)
asc_x  <- posixct_to_decimal(asc_data$datetime, profile_data$datetime)

# Fit predictions
base_slope <- dlmoFit$base_params_fine
base_fitted <- base_slope * (base_x - ip_x) + ip_y
base_residuals <- base_data$melatonin - base_fitted

asc_params <- dlmoFit$ascending_params_fine
if (length(asc_params) == 1) {
  asc_fitted <- asc_params[[1]] * (asc_x - ip_x) + ip_y
} else {
  a <- asc_params$a
  b <- asc_params$b
  c <- asc_params$c
  asc_fitted <- a * asc_x^2 + b * asc_x + c
}
asc_residuals <- asc_data$melatonin - asc_fitted

# --- Residual Bootstrapping ---
n_iter <- 3
dlmo_resid_boot <- numeric(n_iter)

set.seed(42)
for (i in seq_len(n_iter)) {
  print(i)
  res_base <- sample(base_residuals, replace = TRUE)
  res_asc  <- sample(asc_residuals, replace = TRUE)

  boot_base <- base_fitted + res_base
  boot_asc  <- asc_fitted + res_asc

  pseudo_profile <- dplyr::bind_rows(
    dplyr::mutate(base_data, melatonin = boot_base),
    dplyr::mutate(asc_data, melatonin = boot_asc)
  ) |> dplyr::arrange(datetime)

  dlmo_result <- tryCatch({
    calculate_dlmo(data = pseudo_profile)
  }, error = function(e) NULL)

  if (!is.null(dlmo_result)) {
    dlmo_resid_boot[i] <- dlmo_result$ip$inflection_point_fine$x
  } else {
    dlmo_resid_boot[i] <- NA
  }
}

# --- Summary stats ---
dlmo_resid_boot <- dlmo_resid_boot[!is.na(dlmo_resid_boot)]
dlmo_mean <- mean(dlmo_resid_boot)
dlmo_ci <- quantile(dlmo_resid_boot, c(0.025, 0.975))
dlmo_mean_posix <- decimal_to_posixct(dlmo_mean, profile_data$datetime)
dlmo_ci_posix <- decimal_to_posixct(dlmo_ci, profile_data$datetime)
actual_dlmo <- decimal_to_posixct(dlmoFit$inflection_point$x, profile_data$datetime)

# --- Plotting ---
plot_df <- data.frame(dlmo = decimal_to_posixct(dlmo_resid_boot, profile_data$datetime))

density_vals <- density(as.numeric(plot_df$dlmo))
max_hist_count <- max(ggplot_build(
  ggplot(plot_df, aes(x = dlmo)) +
    geom_histogram(binwidth = 60 * 5))$data[[1]]$count)

y_text <- max_hist_count - 0.3
scaled_density <- data.frame(
  x = as.POSIXct(density_vals$x, origin = "1970-01-01", tz = "UTC"),
  y = density_vals$y * max_hist_count / max(density_vals$y)
)

ci_df <- data.frame(
  xmin = dlmo_ci_posix[1],
  xmax = dlmo_ci_posix[2],
  ymin = 0,
  ymax = Inf,
  group = factor("95% CI", levels = c("DLMO estimate", "Bootstrap mean", "95% CI"))
)

line_df <- data.frame(
  label = factor(c("DLMO estimate", "Bootstrap mean"),
                 levels = c("DLMO estimate", "Bootstrap mean", "95% CI")),
  x = c(actual_dlmo, dlmo_mean_posix)
)

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
    title = "Residual Bootstrap: DLMO Estimate Histogram + Scaled Density",
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


# # Convert datetimes to decimal hours (relative to profile_data$datetime)
# decimal_times <- posixct_to_decimal(sample_dlmo$prof$datetime, sample_dlmo$prof$datetime)
#
# # Extract inflection point
# ip_x <- sample_dlmo$inflection_point_fine$x
# ip_y <- sample_dlmo$inflection_point_fine$y
#
# # Subset base and ascending data
# base_data <-  sample_dlmo$prof %>% dplyr::filter(base == 1)
# asc_data <-  sample_dlmo$prof %>% dplyr::filter(ascending == 1)
#
# # Convert their x to decimal
# base_x <- posixct_to_decimal(base_data$datetime, sample_dlmo$datetime)
# asc_x  <- posixct_to_decimal(asc_data$datetime,  sample_dlmo$datetime)
#
# ## --- BASE RESIDUALS (always linear) ---
# base_slope <- sample_dlmo$ip$base_params
# base_predicted <- base_slope * (base_x - ip_x) + ip_y
# base_residuals <- base_data$melatonin - base_predicted
#
# ## --- ASCENDING RESIDUALS (linear or parabolic) ---
# asc_params <- sample_dlmo$ip$ascending_params_fine
#
# if (length(asc_params) == 1) {
#   # Linear case
#   asc_slope <- asc_params[[1]]
#   asc_predicted <- asc_slope * (asc_x - ip_x) + ip_y
# } else {
#   # Parabolic case: y = ax^2 + bx + c
#   a <- asc_params$a
#   b <- asc_params$b
#   c <- asc_params$c
#   asc_predicted <- a * asc_x^2 + b * asc_x + c
# }
#
# asc_residuals <- asc_data$melatonin - asc_predicted
#
# # merge residuals
# combo_residuals <- c(base_residuals, asc_residuals)
#
# # randomly sample residuals with replacement

