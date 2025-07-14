# --- Monte Carlo Bootstrap using Time and Melatonin Noise ---
run_dlmo_mc_bootstrap <- function(sample_dlmo, n_iter = 300, time_sd = NULL, mel_sd = 0.3,
                                  seed = 42, n_workers = parallel::detectCores() - 1) {
  profile_data <- sample_dlmo$prof
  dh_time <- posixct_to_decimal(profile_data$datetime, profile_data$datetime)
  min_interval <- min(diff(sort(dh_time)))
  if (is.null(time_sd)) time_sd <- min_interval / 2

  interp_fun <- splinefun(x = dh_time, y = profile_data$melatonin, method = "monoH.FC")

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
    dlmo_mc <- future_sapply(seq_len(n_iter), future.seed = TRUE, function(i) {
      p()
      library(dlmoR); library(dplyr); library(hms); library(lubridate); library(tibble)

      tryCatch({
        time_pert <- dh_time + rnorm(length(dh_time), mean = 0, sd = time_sd)
        time_sorted <- enforce_min_time_gap(sort(time_pert), min_gap = 1 / 60)

        mel_interp <- interp_fun(time_sorted)
        mel_pert <- mel_interp + rnorm(length(mel_interp), mean = 0, sd = mel_sd)

        df <- tibble(
          datetime = decimal_to_posixct(time_sorted, profile_data$datetime),
          melatonin = mel_pert
        )

        result <- calculate_dlmo(data = df)
        result$ip$inflection_point_fine$x
      }, error = function(e) {
        message(sprintf("❌ MC Iteration %d failed: %s", i, e$message))
        NA_real_
      })
    })
  })

  dlmo_mc <- dlmo_mc[!is.na(dlmo_mc)]
  list(
    bootstrap_values = dlmo_mc,
    mean = mean(dlmo_mc),
    ci = quantile(dlmo_mc, c(0.025, 0.975)),
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
  density_vals <- density(as.numeric(plot_df$dlmo), bw = 60)
  max_hist_count <- max(ggplot_build(
    ggplot(plot_df, aes(x = dlmo)) +
      geom_histogram(binwidth = 60)
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
    geom_histogram(binwidth = 60, fill = "lightblue", color = "black", boundary = 0) +
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

# Example usage:
#boot_mc <- run_dlmo_mc_bootstrap(sample_dlmo, n_iter = 1000)
plot_dlmo_bootstrap(boot_mc, method_label = "Monte Carlo Bootstrap")
