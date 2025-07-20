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

# --- Utilities ---

posixct_to_decimal <- function(posixct_vec, ref = NULL) {
  if (is.null(ref)) ref <- posixct_vec
  as.numeric(difftime(posixct_vec, min(ref), units = "hours"))
}

decimal_to_posixct <- function(decimal_vec, ref_time) {
  origin_time <- min(ref_time)
  origin_time + dhours(decimal_vec)
}

enforce_min_time_gap <- function(time_vec, min_gap = 1 / 60) {
  for (i in 2:length(time_vec)) {
    if ((time_vec[i] - time_vec[i - 1]) < min_gap) {
      time_vec[i] <- time_vec[i - 1] + min_gap
    }
  }
  time_vec
}

# Safely-wrapped noisy DLMO calc
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

    calculate_dlmo(perturbed_df)
  })()
}

safe_calc_dlmo <- safely(calculate_dlmo)

# --- Process single profile ---

process_single_profile <- function(profile_path, k, time_sd, mel_sd, output_dir) {
  profile_name <- tools::file_path_sans_ext(basename(profile_path))
  profile_df <- read_csv(profile_path, show_col_types = FALSE)
  message(glue("📄 {profile_name} — baseline + {k} noisy repetitions"))

  result_list <- list()
  flat_results <- list()

  # 1. Baseline
  baseline_try <- safe_calc_dlmo(profile_df)
  if (!is.null(baseline_try$result)) {
    result <- baseline_try$result
    row <- tibble(
      profile = profile_name,
      type = "baseline",
      rep = NA_integer_,
      dlmo_decimal = result$ip$inflection_point_fine$x,
      dlmo_result = list(result),
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
      time_sd = time_sd,
      mel_sd = mel_sd,
      error_message = baseline_try$error$message
    )
  }
  result_list[[1]] <- list(type = "baseline", result = row)
  flat_results[[1]] <- row

  # 2. Noisy reps
  for (i in seq_len(k)) {
    noise_try <- safe_add_noise_and_calc(profile_df, time_sd, mel_sd)
    if (!is.null(noise_try$result)) {
      result <- noise_try$result
      row <- tibble(
        profile = profile_name,
        type = "noisy",
        rep = i,
        dlmo_decimal = result$ip$inflection_point_fine$x,
        dlmo_result = list(result),
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
        time_sd = time_sd,
        mel_sd = mel_sd,
        error_message = noise_try$error$message
      )
    }
    result_list[[length(result_list) + 1]] <- list(type = "noisy", result = row)
    flat_results[[length(flat_results) + 1]] <- row
  }

  # Save per-profile
  saveRDS(result_list, file = file.path(output_dir, paste0(profile_name, "_dlmo_repeats.rds")))
  flat_results
}

# --- Master pipeline ---

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

  # Save full composite RDS and CSV summary
  saveRDS(all_results_flat, composite_rds)
  write_csv(select(all_results_flat, -dlmo_result), summary_csv)  # CSV excludes bulky objects

  message("✅ All results saved:")
  message("  └─ Per-profile RDS → ", output_dir)
  message("  └─ Composite RDS    → ", composite_rds)
  message("  └─ Summary CSV      → ", summary_csv)
}

process_all_csv_profiles(
  dir_path = "inst/extdata",
  time_sd = 0.005,
  mel_sd = 0.3,
  k = 25,
  output_dir = "noisy_dlmo_results",
  composite_rds = "all_dlmo_results.rds",
  summary_csv = "all_dlmo_summary.csv"
)


# TODO -- plotting!!
