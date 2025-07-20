# DLMO Monte Carlo Noise Simulation Script
# -------------------------------------------
# Author: [Your Name]
# Date: 2025-07-17
# Purpose: Simulate the impact of time sampling jitter and melatonin measurement noise
#          on Dim Light Melatonin Onset (DLMO) estimates across multiple profiles.

# --- Libraries ----------------------------------------------------------
library(readr)      # fast CSV reading
library(dplyr)      # data manipulation
library(purrr)      # functional mapping
library(fs)         # filesystem operations
library(tibble)     # tidy tibbles
library(dlmoR)      # DLMO estimation functions
library(lubridate)  # date-time handling
library(glue)       # string interpolation
library(future)     # parallel plan
library(parallelly) # detect available cores

# --- Utility Functions -------------------------------------------------

posixct_to_decimal <- function(posixct_vec, ref = NULL) {
  if (is.null(ref)) ref <- posixct_vec
  as.numeric(difftime(posixct_vec, min(ref), units = "hours"))
}

decimal_to_posixct <- function(decimal_vec, ref_time) {
  origin_time <- min(ref_time)
  origin_time + dhours(decimal_vec)
}

enforce_min_time_gap <- function(time_vec, min_gap = 1 / 60) {
  for (i in seq(2, length(time_vec))) {
    if ((time_vec[i] - time_vec[i - 1]) < min_gap) {
      time_vec[i] <- time_vec[i - 1] + min_gap
    }
  }
  time_vec
}

# --- Noise and DLMO Wrappers -------------------------------------------

safe_calc_dlmo <- safely(calculate_dlmo)

safe_add_noise_and_calc <- function(profile_df, time_sd, mel_sd) {
  safely(function() {
    dh_time   <- posixct_to_decimal(profile_df$datetime)
    interp_fn <- splinefun(x = dh_time,
                           y = profile_df$melatonin,
                           method = "monoH.FC")

    time_pert   <- dh_time + rnorm(length(dh_time), 0, time_sd)
    time_sorted <- enforce_min_time_gap(sort(time_pert), 1 / 60)

    mel_interp <- interp_fn(time_sorted)
    mel_pert   <- mel_interp + rnorm(length(mel_interp), 0, mel_sd)

    perturbed_df <- tibble(
      datetime   = decimal_to_posixct(time_sorted, profile_df$datetime),
      melatonin  = mel_pert
    )

    list(
      dlmo              = calculate_dlmo(perturbed_df),
      perturbed_profile = perturbed_df
    )
  })()
}

# --- Per-Profile Processor ---------------------------------------------

process_single_profile <- function(profile_path,
                                   k,
                                   time_sd,
                                   mel_sd,
                                   output_dir,
                                   do_baseline = FALSE) {
  profile_name <- tools::file_path_sans_ext(basename(profile_path))
  df           <- read_csv(profile_path, show_col_types = FALSE)
  results      <- list()

  # baseline run (once)
  if (do_baseline) {
    base_try <- run_scenario(df, "baseline", 0, 0)
    results[[1]] <- tibble(
      profile       = profile_name,
      scenario      = "baseline",
      rep           = NA_integer_,
      dlmo_decimal  = if (!is.null(base_try$out$result))
        base_try$out$result$ip$inflection_point_fine$x else NA_real_,
      time_sd       = 0,
      mel_sd        = 0,
      error_message = base_try$out$error$message %||% NA_character_
    )
  }

  # k reps of both_noise and time_only
  idx <- length(results) + 1
  for (i in seq_len(k)) {
    # both_noise
    bn_try <- run_scenario(df, "both_noise", time_sd, mel_sd)
    results[[idx]] <- tibble(
      profile       = profile_name,
      scenario      = "both_noise",
      rep           = i,
      dlmo_decimal  = if (!is.null(bn_try$out$result))
        bn_try$out$result$ip$inflection_point_fine$x else NA_real_,
      time_sd       = time_sd,
      mel_sd        = mel_sd,
      error_message = bn_try$out$error$message %||% NA_character_
    )
    idx <- idx + 1

    # time_only
    to_try <- run_scenario(df, "time_only", time_sd, 0)
    results[[idx]] <- tibble(
      profile       = profile_name,
      scenario      = "time_only",
      rep           = i,
      dlmo_decimal  = if (!is.null(to_try$out$result))
        to_try$out$result$ip$inflection_point_fine$x else NA_real_,
      time_sd       = time_sd,
      mel_sd        = 0,
      error_message = to_try$out$error$message %||% NA_character_
    )
    idx <- idx + 1
  }

  # combine and save
  df_res <- bind_rows(results)
  saveRDS(df_res,
          file = file.path(output_dir, paste0(profile_name, "_results.rds")))
  df_res
}

# helper to dispatch baseline vs noisy
run_scenario <- function(df, scenario, time_sd, mel_sd) {
  if (scenario == "baseline") {
    out <- safe_calc_dlmo(df)
  } else {
    out <- safe_add_noise_and_calc(df, time_sd, mel_sd)
  }
  list(out = out)
}

# --- Main Simulation ---------------------------------------------------

# parameters (edit as needed)
#time_sd_levels <- c(0, 1, 3, 5, 10)  # minutes of jitter to test
time_sd_levels <- c(0, 1, 3)  # minutes of jitter to test
k              <- 2                   # reps per scenario
mel_sd         <- 0.3                 # melatonin noise (ng/mL)

# prepare parallel plan
detectable <- parallelly::availableCores()
n_workers  <- max(1, detectable - 1)
plan(multisession, workers = n_workers)

# iterate over jitter levels
all_profiles <- list.files("inst/extdata", pattern = "\\.csv$",
                           full.names = TRUE)[1:2]

output_base <- "outputs/noisy_dlmo_trials"
dir_create(output_base)

all_results <- list()
for (sd_min in time_sd_levels) {
  t_sd_h <- sd_min / 60
  subdir <- file.path(output_base, paste0(sd_min, "min"))
  dir_create(subdir)

  # baseline only for first level
  do_baseline <- (sd_min == time_sd_levels[1])

  res_df <- future_map_dfr(
    all_profiles,
    process_single_profile,
    k           = k,
    time_sd     = t_sd_h,
    mel_sd      = mel_sd,
    output_dir  = subdir,
    do_baseline = do_baseline
  )

  write_csv(select(res_df, -error_message),
            file.path(subdir, glue("summary_{sd_min}min.csv")))
  all_results[[paste0(sd_min, "min")]] <- res_df
}

# combine all levels into one tibble if desired
combined <- bind_rows(all_results, .id = "time_jitter")
write_csv(combined, file.path(output_base, "all_results_combined.csv"))

message("✅ Simulation complete. Results in ", output_base)
