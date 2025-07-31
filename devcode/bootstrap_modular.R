# dlmo_bootstrap_all.R
# Unified bootstrap framework for DLMO estimation uncertainty in dlmoR

library(dlmoR)
library(dplyr)
library(hms)
library(lubridate)
library(future)
library(future.apply)
library(progressr)
library(ggplot2)


# --- Helper: Enforce minimum gap for sorted times ---
enforce_min_time_gap <- function(time_vec, min_gap = 1 / 60) {
  for (i in 2:length(time_vec)) {
    if ((time_vec[i] - time_vec[i - 1]) < min_gap) {
      time_vec[i] <- time_vec[i - 1] + min_gap
    }
  }
  time_vec
}

# --- Unified Bootstrap Dispatcher ---
dlmo_bootstrap <- function(sample_dlmo, method = c("monte_carlo", "residual", "wild", "hybrid"),
                           n_iter = 300, time_sd = NULL, mel_cv = 0.079,
                           wild_type = c("rademacher", "normal"),
                           clip_negatives = TRUE, seed = 42,
                           n_workers = parallel::detectCores() - 1) {

  method <- match.arg(method)
  wild_type <- match.arg(wild_type)

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

  # Parallel setup
  plan(multisession, workers = n_workers)
  handlers(global = TRUE)
  handlers("rstudio")
  set.seed(seed)

  with_progress({
    p <- progressor(steps = n_iter)
    results <- future_sapply(seq_len(n_iter), future.seed = TRUE, function(i) {
      p()
      library(dlmoR); library(dplyr); library(hms); library(lubridate)

      tryCatch({
        if (method == "monte_carlo") {
          jittered_base <- base_x + rnorm(length(base_x), 0, time_sd)
          jittered_asc  <- asc_x  + rnorm(length(asc_x),  0, time_sd)
          jittered_base <- enforce_min_time_gap(sort(jittered_base))
          jittered_asc  <- enforce_min_time_gap(sort(jittered_asc))

          mel_base <- base_data$melatonin * rnorm(length(base_x), 1, mel_cv)
          mel_asc  <- asc_data$melatonin * rnorm(length(asc_x),  1, mel_cv)

          if (clip_negatives) {
            mel_base[mel_base < 0] <- 0
            mel_asc[mel_asc < 0] <- 0
          }

          pseudo_profile <- bind_rows(
            tibble(datetime = decimal_to_posixct(jittered_base, datetime_ref), melatonin = mel_base),
            tibble(datetime = decimal_to_posixct(jittered_asc, datetime_ref),  melatonin = mel_asc)
          ) |> arrange(datetime)

        } else if (method == "residual") {
          res_base <- sample(base_residuals, replace = TRUE)
          res_asc  <- sample(asc_residuals, replace = TRUE)

          mel_base <- base_fitted + res_base
          mel_asc  <- asc_fitted + res_asc

          pseudo_profile <- bind_rows(
            mutate(base_data, melatonin = mel_base),
            mutate(asc_data,  melatonin = mel_asc)
          ) |> arrange(datetime)

        } else if (method == "wild") {
          w_base <- if (wild_type == "rademacher") sample(c(-1, 1), length(base_residuals), TRUE) else rnorm(length(base_residuals))
          w_asc  <- if (wild_type == "rademacher") sample(c(-1, 1), length(asc_residuals), TRUE) else rnorm(length(asc_residuals))

          mel_base <- base_fitted + base_residuals * w_base
          mel_asc  <- asc_fitted + asc_residuals * w_asc

          pseudo_profile <- bind_rows(
            mutate(base_data, melatonin = mel_base),
            mutate(asc_data,  melatonin = mel_asc)
          ) |> arrange(datetime)

        } else if (method == "hybrid") {
          jittered_base <- enforce_min_time_gap(sort(base_x + rnorm(length(base_x), 0, time_sd)))
          jittered_asc  <- enforce_min_time_gap(sort(asc_x  + rnorm(length(asc_x),  0, time_sd)))

          res_base <- sample(base_residuals, replace = TRUE)
          res_asc  <- sample(asc_residuals, replace = TRUE)

          mel_base <- base_slope * (jittered_base - ip_x) + ip_y + res_base
          mel_asc <- if (length(asc_params) == 1) {
            asc_params[[1]] * (jittered_asc - ip_x) + ip_y + res_asc
          } else {
            a <- asc_params$a; b <- asc_params$b; c <- asc_params$c
            a * jittered_asc^2 + b * jittered_asc + c + res_asc
          }

          pseudo_profile <- bind_rows(
            tibble(datetime = decimal_to_posixct(jittered_base, datetime_ref), melatonin = mel_base),
            tibble(datetime = decimal_to_posixct(jittered_asc,  datetime_ref), melatonin = mel_asc)
          ) |> arrange(datetime)
        }

        pseudo_profile$time <- as_hms(as_datetime(pseudo_profile$datetime))
        result <- calculate_dlmo(data = pseudo_profile)
        result$ip$inflection_point_fine$x

      }, error = function(e) {
        message(sprintf("❌ Iteration %d failed: %s", i, e$message))
        NA_real_
      })
    })
  })

  results <- results[!is.na(results)]
  list(
    bootstrap_values = results,
    mean = mean(results),
    ci = quantile(results, c(0.025, 0.975)),
    actual = sample_dlmo$ip$inflection_point_fine$x,
    ref_time = datetime_ref,
    dlmo_label = sample_dlmo$dlmo$fine$time
  )
}

# --- Plotting Function ---
plot_dlmo_bootstrap <- function(boot_result, method_label = "Bootstrap", bw = 300) {
  dlmo_resid_boot <- boot_result$bootstrap_values
  dlmo_mean <- boot_result$mean
  dlmo_ci <- boot_result$ci
  actual_dlmo <- decimal_to_posixct(boot_result$actual, boot_result$ref_time)
  dlmo_mean_posix <- decimal_to_posixct(dlmo_mean, boot_result$ref_time)
  dlmo_ci_posix <- decimal_to_posixct(dlmo_ci, boot_result$ref_time)
  format(dlmo_ci_posix, "%H:%M:%S")

  plot_df <- data.frame(dlmo = decimal_to_posixct(dlmo_resid_boot, boot_result$ref_time))
  density_vals <- density(as.numeric(plot_df$dlmo), bw = bw)
  max_hist_count <- max(ggplot_build(ggplot(plot_df, aes(x = dlmo)) +
                                       geom_histogram(binwidth = bw))$data[[1]]$count)

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

  ggplot(plot_df, aes(x = dlmo)) +
    geom_rect(data = ci_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = group),
              alpha = 0.2, inherit.aes = FALSE) +
    geom_histogram(binwidth = bw, fill = "darkgrey", color = "black", boundary = 0) +
    geom_line(data = scaled_density, aes(x = x, y = y), color = "darkslategrey", linewidth = 1) +
    geom_vline(data = line_df, aes(xintercept = x, color = label, linetype = label), linewidth = 1) +
    scale_color_manual(values = c("DLMO estimate" = "deeppink3", "Bootstrap mean" = "darkgoldenrod")) +
    scale_linetype_manual(values = c("DLMO estimate" = "solid", "Bootstrap mean" = "solid")) +
    scale_fill_manual(name = "Reference", values = c("95% CI" = "cadetblue3")) +
    scale_x_datetime(date_labels = "%H:%M", date_breaks = "15 min", timezone = "UTC", expand = c(0.01,0.01)) +
    scale_y_continuous(
      name = "Count"
      #sec.axis = sec_axis(~ . * max(density_vals$y) / max_hist_count, name = "Density")
    ) +
    labs(
      title = paste0("DLMO estimate: ", method_label),
      x = "Time [hh:mm]"
    ) +
    # annotate("text", x = actual_dlmo, y = y_text,
    #          label = paste0("DLMO: ", boot_result$dlmo_label),
    #          color = "deeppink", hjust = -0.1, vjust = 1, size = 3) +
    # annotate("text", x = dlmo_mean_posix, y = y_text,
    #          label = paste0("Mean: ", format(dlmo_mean_posix, "%H:%M:%S")),
    #          color = "blue", hjust = -0.1, vjust = 2.5, size = 3) +
    annotation_custom(
      grid::textGrob(
        label = paste0("DLMO: ", boot_result$dlmo_label),
        x = unit(0.02, "npc"), y = unit(0.98, "npc"),
        just = c("left", "top"),
        gp = gpar(col = "deeppink3", fontsize = 10, fontface = "bold")
      )
    ) +
    annotation_custom(
      grid::textGrob(
        label = paste0("Mean: ", format(dlmo_mean_posix, "%H:%M:%S")),
        x = unit(0.02, "npc"), y = unit(0.93, "npc"),
        just = c("left", "top"),
        gp = gpar(col = "darkgoldenrod", fontsize = 10, fontface = "bold")
      )
    ) +
    annotation_custom(
      grid::textGrob(
        label = paste0("95% CI: ", paste(format(dlmo_ci_posix, "%H:%M:%S"), collapse = " – ")),
        x = unit(0.02, "npc"), y = unit(0.88, "npc"),
        just = c("left", "top"),
        gp = gpar(col = "cadetblue4", fontsize = 10, fontface = "bold")
      )
    )+
    theme_minimal() +
    guides(color = guide_legend(order = 1),
           linetype = guide_legend(order = 1),
           fill = guide_legend(order = 2)) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_blank(),
      legend.box = "horizontal",
      legend.margin = margin(t = 0),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )
}
