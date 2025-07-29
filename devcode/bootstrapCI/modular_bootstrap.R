run_dlmo_bootstrap <- function(sample_dlmo, method = c("montecarlo", "residual", "wild", "hybrid"),
                               n_iter = 300,
                               time_sd = NULL,
                               mel_cv = 0.079,
                               wild_type = c("rademacher", "normal"),
                               clip_negatives = TRUE,
                               seed = 42,
                               n_workers = parallel::detectCores() - 1,
                               plot = FALSE) {
  method <- match.arg(method)
  wild_type <- match.arg(wild_type)

  stopifnot("ip" %in% names(sample_dlmo))
  stopifnot("prof" %in% names(sample_dlmo))

  library(future)
  library(future.apply)
  library(progressr)
  library(dplyr)
  library(hms)
  library(lubridate)

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

  # for Monte Carlo: estimate default time_sd
  if (is.null(time_sd)) {
    min_interval <- min(diff(sort(decimal_times)))
    time_sd <- min_interval / 2
  }

  plan(multisession, workers = n_workers)
  handlers(global = TRUE)
  handlers("rstudio")
  set.seed(seed)

  enforce_min_time_gap <- function(time_vec, min_gap = 1 / 60) {
    for (i in 2:length(time_vec)) {
      gap <- time_vec[i] - time_vec[i - 1]
      if (gap < min_gap) {
        time_vec[i] <- time_vec[i - 1] + min_gap
      }
    }
    time_vec
  }

  with_progress({
    p <- progressor(steps = n_iter)
    dlmo_boot <- future_sapply(seq_len(n_iter), future.seed = TRUE, function(i) {
      p()
      library(dlmoR); library(dplyr); library(hms); library(lubridate)

      tryCatch({
        if (method == "residual") {
          res_base <- sample(base_residuals, replace = TRUE)
          res_asc  <- sample(asc_residuals,  replace = TRUE)
          boot_base <- base_fitted + res_base
          boot_asc  <- asc_fitted + res_asc
          new_times_base <- base_data$datetime
          new_times_asc <- asc_data$datetime

        } else if (method == "wild") {
          w_base <- if (wild_type == "rademacher") sample(c(-1, 1), length(base_residuals), replace = TRUE) else rnorm(length(base_residuals))
          w_asc  <- if (wild_type == "rademacher") sample(c(-1, 1), length(asc_residuals),  replace = TRUE) else rnorm(length(asc_residuals))
          boot_base <- base_fitted + base_residuals * w_base
          boot_asc  <- asc_fitted + asc_residuals * w_asc
          new_times_base <- base_data$datetime
          new_times_asc <- asc_data$datetime

        } else if (method == "montecarlo") {
          base_mel <- base_data$melatonin * (1 + rnorm(length(base_data$melatonin), mean = 0, sd = mel_cv))
          asc_mel  <- asc_data$melatonin * (1 + rnorm(length(asc_data$melatonin), mean = 0, sd = mel_cv))

          base_jittered <- enforce_min_time_gap(sort(base_x + rnorm(length(base_x), 0, time_sd)))
          asc_jittered  <- enforce_min_time_gap(sort(asc_x + rnorm(length(asc_x),  0, time_sd)))

          new_times_base <- decimal_to_posixct(base_jittered, datetime_ref)
          new_times_asc  <- decimal_to_posixct(asc_jittered, datetime_ref)

          boot_base <- base_mel
          boot_asc  <- asc_mel

        } else if (method == "hybrid") {
          base_jittered <- enforce_min_time_gap(sort(base_x + rnorm(length(base_x), 0, time_sd)))
          asc_jittered  <- enforce_min_time_gap(sort(asc_x + rnorm(length(asc_x),  0, time_sd)))
          res_base <- sample(base_residuals, replace = TRUE)
          res_asc  <- sample(asc_residuals,  replace = TRUE)
          boot_base <- base_slope * (base_jittered - ip_x) + ip_y + res_base
          if (length(asc_params) == 1) {
            boot_asc <- asc_params[[1]] * (asc_jittered - ip_x) + ip_y + res_asc
          } else {
            a <- asc_params$a; b <- asc_params$b; c <- asc_params$c
            boot_asc <- a * asc_jittered^2 + b * asc_jittered + c + res_asc
          }
          new_times_base <- decimal_to_posixct(base_jittered, datetime_ref)
          new_times_asc  <- decimal_to_posixct(asc_jittered,  datetime_ref)
        }

        pseudo_profile <- bind_rows(
          mutate(base_data, datetime = new_times_base, melatonin = boot_base),
          mutate(asc_data,  datetime = new_times_asc,  melatonin = boot_asc)
        ) |> arrange(datetime)

        if (clip_negatives) {
          pseudo_profile$melatonin[pseudo_profile$melatonin < 0] <- 0
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

  dlmo_boot <- dlmo_boot[!is.na(dlmo_boot)]
  result <- list(
    bootstrap_values = dlmo_boot,
    mean = mean(dlmo_boot),
    ci = quantile(dlmo_boot, c(0.025, 0.975)),
    actual = sample_dlmo$ip$inflection_point_fine$x,
    ref_time = datetime_ref,
    dlmo_label = sample_dlmo$dlmo$fine$time,
    method = method
  )

  if (plot) {
    plot_dlmo_bootstrap(result)
  }

  return(result)
}
