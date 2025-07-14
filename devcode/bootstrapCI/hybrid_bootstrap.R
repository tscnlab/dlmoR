# --- Hybrid Bootstrap: Jittered Time + Resampled Residuals ---
run_dlmo_hybrid_bootstrap <- function(sample_dlmo, n_iter = 300, time_sd = NULL,
                                      seed = 42, n_workers = parallel::detectCores() - 1) {
  stopifnot("ip" %in% names(sample_dlmo))
  stopifnot("prof" %in% names(sample_dlmo))

  profile_data <- sample_dlmo$prof
  ip_x <- sample_dlmo$ip$inflection_point_fine$x
  ip_y <- sample_dlmo$ip$inflection_point_fine$y
  datetime_ref <- profile_data$datetime
  decimal_times <- posixct_to_decimal(datetime_ref, datetime_ref)

  min_interval <- min(diff(sort(decimal_times)))
  if (is.null(time_sd)) time_sd <- min_interval / 2

  first_base_idx <- which(profile_data$base == 1)[1]
  last_asc_idx <- tail(which(profile_data$ascending == 1), 1)
  middle_range <- profile_data[first_base_idx:last_asc_idx, ]
  middle_times <- decimal_times[first_base_idx:last_asc_idx]
  base_data <- middle_range[middle_times <= ip_x, ]
  asc_data  <- middle_range[middle_times >= ip_x, ]

  base_x <- posixct_to_decimal(base_data$datetime, datetime_ref)
  asc_x  <- posixct_to_decimal(asc_data$datetime, datetime_ref)

  base_slope <- sample_dlmo$ip$base_params_fine
  base_fitted <- base_slope * (base_x - ip_x) + ip_y
  base_residuals <- base_data$melatonin - base_fitted

  asc_params <- sample_dlmo$ip$ascending_params_fine
  asc_fitted <- if (length(asc_params) == 1) {
    asc_params[[1]] * (asc_x - ip_x) + ip_y
  } else {
    a <- asc_params$a; b <- asc_params$b; c <- asc_params$c
    a * asc_x^2 + b * asc_x + c
  }
  asc_residuals <- asc_data$melatonin - asc_fitted

  enforce_min_time_gap <- function(time_vec, min_gap = 1 / 60) {
    for (i in 2:length(time_vec)) {
      gap <- time_vec[i] - time_vec[i - 1]
      if (gap < min_gap) {
        time_vec[i] <- time_vec[i - 1] + min_gap
      }
    }
    time_vec
  }

  set.seed(seed)
  plan(multisession, workers = n_workers)
  handlers(global = TRUE)
  handlers("rstudio")

  with_progress({
    p <- progressor(steps = n_iter)
    dlmo_hybrid <- future_sapply(seq_len(n_iter), future.seed = TRUE, function(i) {
      p()
      library(dlmoR); library(dplyr); library(hms); library(lubridate)

      tryCatch({
        # Jitter time
        base_jittered <- base_x + rnorm(length(base_x), 0, time_sd)
        asc_jittered  <- asc_x +  rnorm(length(asc_x),  0, time_sd)

        # Enforce time monotonicity
        base_jittered <- enforce_min_time_gap(sort(base_jittered))
        asc_jittered  <- enforce_min_time_gap(sort(asc_jittered))

        # Resample residuals
        res_base <- sample(base_residuals, replace = TRUE)
        res_asc  <- sample(asc_residuals,  replace = TRUE)

        # Reconstruct signal
        base_new <- base_slope * (base_jittered - ip_x) + ip_y + res_base
        asc_new  <- if (length(asc_params) == 1) {
          asc_params[[1]] * (asc_jittered - ip_x) + ip_y + res_asc
        } else {
          a <- asc_params$a; b <- asc_params$b; c <- asc_params$c
          a * asc_jittered^2 + b * asc_jittered + c + res_asc
        }

        pseudo_profile <- bind_rows(
          tibble(datetime = decimal_to_posixct(base_jittered, datetime_ref), melatonin = base_new),
          tibble(datetime = decimal_to_posixct(asc_jittered,  datetime_ref), melatonin = asc_new)
        ) |> arrange(datetime)

        pseudo_profile$time <- as_hms(as_datetime(pseudo_profile$datetime))
        result <- calculate_dlmo(data = pseudo_profile)
        result$ip$inflection_point_fine$x
      }, error = function(e) {
        message(sprintf("❌ Iteration %d failed: %s", i, e$message))
        NA_real_
      })
    })
  })

  dlmo_hybrid <- dlmo_hybrid[!is.na(dlmo_hybrid)]
  list(
    bootstrap_values = dlmo_hybrid,
    mean = mean(dlmo_hybrid),
    ci = quantile(dlmo_hybrid, c(0.025, 0.975)),
    actual = sample_dlmo$ip$inflection_point_fine$x,
    ref_time = profile_data$datetime,
    dlmo_label = sample_dlmo$dlmo$fine$time
  )
}

plot_dlmo_bootstrap <- function(boot_result, method_label = "Bootstrap") {
  # Extract components
  dlmo_resid_boot <- boot_result$bootstrap_values
  dlmo_mean <- boot_result$mean
  dlmo_ci <- boot_result$ci
  actual_dlmo <- decimal_to_posixct(boot_result$actual, boot_result$ref_time)
  dlmo_mean_posix <- decimal_to_posixct(dlmo_mean, boot_result$ref_time)
  dlmo_ci_posix <- decimal_to_posixct(dlmo_ci, boot_result$ref_time)

  # Data for plot
  plot_df <- data.frame(dlmo = decimal_to_posixct(dlmo_resid_boot, boot_result$ref_time))
  density_vals <- density(as.numeric(plot_df$dlmo))
  max_hist_count <- max(ggplot_build(
    ggplot(plot_df, aes(x = dlmo)) +
      geom_histogram(binwidth = 60*5)
  )$data[[1]]$count)

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

  y_text <- max_hist_count - 0.3

  # Plot
  ggplot(plot_df, aes(x = dlmo)) +
    geom_histogram(binwidth = 60*5, fill = "lightblue", color = "black", boundary = 0) +
    geom_rect(data = ci_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = group),
              alpha = 0.15, inherit.aes = FALSE) +
    geom_line(data = scaled_density, aes(x = x, y = y), color = "darkblue", linewidth = 1) +
    geom_vline(data = line_df, aes(xintercept = x, color = label, linetype = label), linewidth = 1) +
    scale_color_manual(values = c("DLMO estimate" = "deeppink", "Bootstrap mean" = "blue")) +
    scale_linetype_manual(values = c("DLMO estimate" = "dotdash", "Bootstrap mean" = "solid")) +
    scale_fill_manual(name = "Reference", values = c("95% CI" = "azure4")) +
    scale_x_datetime(date_labels = "%H:%M", date_breaks = "15 min", timezone = "UTC", expand = c(0.01, 0.01)) +
    scale_y_continuous(
      name = "Count",
      sec.axis = sec_axis(~ . * max(density_vals$y) / max_hist_count, name = "Density")
    ) +
    labs(
      title = paste0(method_label, ": DLMO Estimate Histogram + Scaled Density"),
      x = "DLMO Estimate [hh:mm]"
    ) +
    annotate("text", x = actual_dlmo, y = y_text,
             label = paste0("DLMO: ", boot_result$dlmo_label),
             color = "deeppink", hjust = -0.1, vjust = 1, size = 3) +
    annotate("text", x = dlmo_mean_posix, y = y_text,
             label = paste0("Mean: ", format(dlmo_mean_posix, "%H:%M:%S")),
             color = "blue", hjust = -0.1, vjust = 2.5, size = 3) +
    theme_minimal() +
    guides(color = guide_legend(order = 1),
           linetype = guide_legend(order = 1),
           fill = guide_legend(order = 2)) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_blank(),
      legend.box = "horizontal",
      legend.margin = margin(t = 0)
    )
}

# --- Example: Run and Plot ---
# filename <- system.file("extdata/sample_melatonin_profile.csv", package = "dlmoR")
# sample_dlmo <- calculate_dlmo(file_path = filename, threshold = 2.3)
# boot_hybrid <- run_dlmo_hybrid_bootstrap(sample_dlmo, n_iter = 30)
plot_dlmo_bootstrap(boot_hybrid, method_label = "Hybrid Bootstrap")

# --- Validation Plot: Original, Interpolated, Hybrid Sample ---
dh_time <- posixct_to_decimal(sample_dlmo$prof$datetime, sample_dlmo$prof$datetime)
interp_fun <- splinefun(x = dh_time, y = sample_dlmo$prof$melatonin, method = "monoH.FC")
time_grid <- seq(min(dh_time), max(dh_time), length.out = 500)
mel_interp_grid <- interp_fun(time_grid)
time_sorted <- sort(dh_time + rnorm(length(dh_time), mean = 0, sd = min(diff(dh_time)) / 2))
interp_vals <- interp_fun(time_sorted)

# Residual resampling
ip_x <- sample_dlmo$ip$inflection_point_fine$x
ip_y <- sample_dlmo$ip$inflection_point_fine$y
profile_data <- sample_dlmo$prof
base_idx <- which(profile_data$base == 1)[1]
asc_idx <- tail(which(profile_data$ascending == 1), 1)
base_data <- profile_data[base_idx:asc_idx, ] |> filter(posixct_to_decimal(datetime, datetime) <= ip_x)
asc_data <- profile_data[base_idx:asc_idx, ] |> filter(posixct_to_decimal(datetime, datetime) >= ip_x)

base_x <- posixct_to_decimal(base_data$datetime, sample_dlmo$prof$datetime)
asc_x <- posixct_to_decimal(asc_data$datetime, sample_dlmo$prof$datetime)
base_slope <- sample_dlmo$ip$base_params_fine
base_fit <- base_slope * (base_x - ip_x) + ip_y
base_res <- base_data$melatonin - base_fit

asc_params <- sample_dlmo$ip$ascending_params_fine
asc_fit <- if (length(asc_params) == 1) {
  asc_params[[1]] * (asc_x - ip_x) + ip_y
} else {
  a <- asc_params$a; b <- asc_params$b; c <- asc_params$c
  a * asc_x^2 + b * asc_x + c
}
asc_res <- asc_data$melatonin - asc_fit

# One hybrid sample
res_base <- sample(base_res, replace = TRUE)
res_asc  <- sample(asc_res,  replace = TRUE)
mel_base <- base_fit + res_base
mel_asc  <- asc_fit + res_asc

# Jittered time
time_base <- sort(base_x + rnorm(length(base_x), mean = 0, sd = min(diff(dh_time)) / 2))
time_asc  <- sort(asc_x + rnorm(length(asc_x),  mean = 0, sd = min(diff(dh_time)) / 2))

hybrid_df <- bind_rows(
  tibble(datetime = decimal_to_posixct(time_base, sample_dlmo$prof$datetime), melatonin = mel_base, source = "Hybrid sample"),
  tibble(datetime = decimal_to_posixct(time_asc, sample_dlmo$prof$datetime), melatonin = mel_asc, source = "Hybrid sample")
)

# Plot
original_df <- sample_dlmo$prof |> mutate(source = "Original measurements")
interp_df <- tibble(datetime = decimal_to_posixct(time_grid, sample_dlmo$prof$datetime), melatonin = mel_interp_grid, source = "Interpolated curve")

combined_df <- bind_rows(original_df, interp_df, hybrid_df)

library(ggplot2)
ggplot(combined_df, aes(x = datetime, y = melatonin, color = source)) +
  geom_line(data = filter(combined_df, source == "Interpolated curve"), linewidth = 1) +
  geom_point(data = filter(combined_df, source != "Interpolated curve"), size = 2, alpha = 0.8) +
  scale_color_manual(
    values = c(
      "Original measurements" = "black",
      "Interpolated curve" = "blue",
      "Hybrid sample" = "darkred"
    )
  ) +
  labs(
    title = "Validation of Hybrid Bootstrap Sample",
    x = "Time",
    y = "Melatonin (pg/mL)",
    color = "Data Source"
  ) +
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "1 hour") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )
