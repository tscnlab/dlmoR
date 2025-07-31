run_dlmo_wild_bootstrap <- function(sample_dlmo, n_iter = 300, wild_type = c("rademacher", "normal"),
                                    seed = 42, n_workers = parallel::detectCores() - 1) {
  wild_type <- match.arg(wild_type)

  profile_data <- sample_dlmo$prof
  ip_x <- sample_dlmo$ip$inflection_point_fine$x
  ip_y <- sample_dlmo$ip$inflection_point_fine$y
  datetime_ref <- profile_data$datetime
  decimal_times <- posixct_to_decimal(datetime_ref, datetime_ref)

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

  plan(multisession, workers = n_workers)
  handlers(global = TRUE)
  handlers("rstudio")
  set.seed(seed)

  with_progress({
    p <- progressor(steps = n_iter)
    dlmo_boot <- future_sapply(seq_len(n_iter), future.seed = TRUE, function(i) {
      p()
      library(dlmoR); library(dplyr); library(hms); library(lubridate)

      tryCatch({
        # Wild weights
        if (wild_type == "rademacher") {
          w_base <- sample(c(-1, 1), length(base_residuals), replace = TRUE)
          w_asc  <- sample(c(-1, 1), length(asc_residuals), replace = TRUE)
        } else if (wild_type == "normal") {
          w_base <- rnorm(length(base_residuals), mean = 0, sd = 1)
          w_asc  <- rnorm(length(asc_residuals), mean = 0, sd = 1)
        }

        boot_base <- base_fitted + base_residuals * w_base
        boot_asc  <- asc_fitted + asc_residuals * w_asc

        pseudo_profile <- bind_rows(
          mutate(base_data, melatonin = boot_base),
          mutate(asc_data, melatonin = boot_asc)
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

  dlmo_boot <- dlmo_boot[!is.na(dlmo_boot)]
  list(
    bootstrap_values = dlmo_boot,
    mean = mean(dlmo_boot),
    ci = quantile(dlmo_boot, c(0.025, 0.975)),
    actual = sample_dlmo$ip$inflection_point_fine$x,
    ref_time = datetime_ref,
    dlmo_label = sample_dlmo$dlmo$fine$time
  )
}

filename <- system.file("extdata/civibe_melatonin_FD207_day2.csv", package = "dlmoR")
sample_dlmo <- calculate_dlmo(file_path = filename, threshold = 2.3)

boot_rademacher <- run_dlmo_wild_bootstrap(sample_dlmo, n_iter = 300, wild_type = "rademacher")
#boot_normal     <- run_dlmo_wild_bootstrap(sample_dlmo, n_iter = 300, wild_type = "normal")


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
  density_vals <- density(as.numeric(plot_df$dlmo), bw=60)
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

plot_dlmo_bootstrap(boot_rademacher, method_label = "Wild Bootstrap (Rademacher)")
plot_dlmo_bootstrap(boot_normal, method_label = "Wild Bootstrap (Normal)")
