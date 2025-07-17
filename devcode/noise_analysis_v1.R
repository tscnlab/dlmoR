# DLMO Simulation Pipeline with Time and Melatonin Noise
# This pipeline processes melatonin profiles to estimate DLMO using both a baseline
# and multiple noisy perturbations. It simulates irregular sampling by adding time jitter,
# and adds Gaussian noise to melatonin concentration. Perturbed profiles are saved
# alongside results for traceability.

library(readr)
library(dplyr)
library(purrr)
library(fs)
library(tibble)
library(dlmoR)
library(lubridate)
library(glue)
library(future)
library(furrr)

# --- Utility Functions ---

# Convert POSIXct datetime to decimal hours since start
posixct_to_decimal <- function(posixct_vec, ref = NULL) {
  if (is.null(ref)) ref <- posixct_vec
  as.numeric(difftime(posixct_vec, min(ref), units = "hours"))
}

# Convert decimal hours back to POSIXct datetime
decimal_to_posixct <- function(decimal_vec, ref_time) {
  origin_time <- min(ref_time)
  origin_time + dhours(decimal_vec)
}

# Enforce a minimum time gap between successive time points (e.g. 1 minute)
enforce_min_time_gap <- function(time_vec, min_gap = 1 / 60) {
  for (i in 2:length(time_vec)) {
    if ((time_vec[i] - time_vec[i - 1]) < min_gap) {
      time_vec[i] <- time_vec[i - 1] + min_gap
    }
  }
  time_vec
}

# --- DLMO Estimation with Noise ---

# Apply Gaussian jitter to time, noise to melatonin, and calculate DLMO
# Returns both DLMO result and perturbed profile
safe_add_noise_and_calc <- function(profile_df, time_sd, mel_sd) {
  safely(function() {
    dh_time <- posixct_to_decimal(profile_df$datetime)
    interp_fun <- splinefun(x = dh_time, y = profile_df$melatonin, method = "monoH.FC")

    time_pert <- dh_time + rnorm(length(dh_time), 0, time_sd)
    time_sorted <- enforce_min_time_gap(sort(time_pert), min_gap = 1 / 60)

    mel_interp <- interp_fun(time_sorted)
    mel_pert <- mel_interp + rnorm(length(mel_interp), 0, mel_sd)

    perturbed_df <- tibble(
      datetime = decimal_to_posixct(time_sorted, profile_df$datetime),
      melatonin = mel_pert
    )

    list(
      dlmo = calculate_dlmo(perturbed_df),
      perturbed_profile = perturbed_df
    )
  })()
}

# Safe wrapper around baseline DLMO estimation
safe_calc_dlmo <- safely(calculate_dlmo)

# --- Process Single Profile ---

# Runs baseline + k noisy DLMO repetitions on a single profile file
# Saves all results (including perturbed profiles) as a .rds file
process_single_profile <- function(profile_path, k, time_sd, mel_sd, output_dir) {
  profile_name <- tools::file_path_sans_ext(basename(profile_path))
  profile_df <- read_csv(profile_path, show_col_types = FALSE)
  message(glue("\U0001F4C4 {profile_name} — baseline + {k} noisy repetitions"))

  result_list <- list()
  flat_results <- list()

  # 1. Baseline (no noise)
  baseline_try <- safe_calc_dlmo(profile_df)
  if (!is.null(baseline_try$result)) {
    result <- baseline_try$result
    row <- tibble(
      profile = profile_name,
      type = "baseline",
      rep = NA_integer_,
      dlmo_decimal = result$ip$inflection_point_fine$x,
      dlmo_result = list(result),
      perturbed_profile = list(NULL),  # no perturbation for baseline
      time_sd = time_sd,
      mel_sd = mel_sd,
      error_message = NA_character_
    )
  } else {
    row <- tibble(
      profile = profile_name,
      type = "baseline",
      rep = NA_integer_,
      dlmo_decimal = NA_real_,
      dlmo_result = list(NULL),
      perturbed_profile = list(NULL),
      time_sd = time_sd,
      mel_sd = mel_sd,
      error_message = baseline_try$error$message
    )
  }
  result_list[[1]] <- list(type = "baseline", result = row)
  flat_results[[1]] <- row

  # 2. Noisy repetitions
  for (i in seq_len(k)) {
    noise_try <- safe_add_noise_and_calc(profile_df, time_sd, mel_sd)
    if (!is.null(noise_try$result$dlmo)) {
      result <- noise_try$result$dlmo
      perturbed_df <- noise_try$result$perturbed_profile
      row <- tibble(
        profile = profile_name,
        type = "noisy",
        rep = i,
        dlmo_decimal = result$ip$inflection_point_fine$x,
        dlmo_result = list(result),
        perturbed_profile = list(perturbed_df),
        time_sd = time_sd,
        mel_sd = mel_sd,
        error_message = NA_character_
      )
    } else {
      row <- tibble(
        profile = profile_name,
        type = "noisy",
        rep = i,
        dlmo_decimal = NA_real_,
        dlmo_result = list(NULL),
        perturbed_profile = list(NULL),
        time_sd = time_sd,
        mel_sd = mel_sd,
        error_message = noise_try$error$message
      )
    }
    result_list[[length(result_list) + 1]] <- list(type = "noisy", result = row)
    flat_results[[length(flat_results) + 1]] <- row
  }

  # Save per-profile results
  saveRDS(result_list, file = file.path(output_dir, paste0(profile_name, "_dlmo_repeats.rds")))
  flat_results
}

# --- Process All Profiles in Directory ---

# Applies the above processing across all CSV files in a directory
# Saves a composite .rds and a summary .csv
process_all_csv_profiles <- function(
    dir_path = "inst/extdata",
    time_sd = 0.005,
    mel_sd = 0.3,
    k = 25,
    output_dir = "noisy_dlmo_results",
    composite_rds = "all_dlmo_results.rds",
    summary_csv = "all_dlmo_summary.csv"
) {
  dir_create(output_dir)
  profile_files <- list.files(dir_path, pattern = "\\.csv$", full.names = TRUE)

  plan(multisession, workers = parallel::detectCores() - 1)

  all_results_flat <- future_map(
    profile_files,
    ~process_single_profile(.x, k, time_sd, mel_sd, output_dir)
  ) %>%
    flatten_df()

  saveRDS(all_results_flat, composite_rds)
  write_csv(select(all_results_flat, -dlmo_result, -perturbed_profile), summary_csv)

  message("\u2705 All results saved:")
  message("  \u2514\u2500 Per-profile RDS \u2192 ", output_dir)
  message("  \u2514\u2500 Composite RDS    \u2192 ", composite_rds)
  message("  \u2514\u2500 Summary CSV      \u2192 ", summary_csv)
}
