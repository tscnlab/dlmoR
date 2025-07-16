# ------------------------------------------------------------------------------
# DLMO Resampling Sensitivity Analysis
#
# This script evaluates how estimated DLMO (Dim Light Melatonin Onset) values
# change as a function of temporal sampling interval. It loads multiple melatonin
# profiles, computes the original DLMO for each, then resamples each profile at
# various intervals using PCHIP interpolation. DLMO is re-estimated for each
# resampled profile, and the differences are recorded and visualized.
#
# Key features:
# - Parallel processing across profiles and resample intervals (via {furrr})
# - Error handling at both profile and interval levels
# - Progress bar support via {progressr}
# - Output includes results table, error logs, and summary plots
#
# Outputs:
# - RDS file of all successful results
# - CSV file of all errors
# - Two plots:
#     1. DLMO vs. Sampling Interval
#     2. Change in DLMO vs. Sampling Interval
# ------------------------------------------------------------------------------

# -----------------------------
# Required packages
# -----------------------------
library(readr)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pracma)    # for pchip interpolation
library(tools)     # for file name handling
library(furrr)     # for parallelization
library(future)    # for parallelization
library(progressr) # for progress bar

# -----------------------------
# 1a. Load all melatonin profiles from CSVs
# -----------------------------
profile_folder <- system.file("extdata/", package = "dlmoR")
csv_files <- list.files(profile_folder, pattern = "\\.csv$", full.names = TRUE)

profiles <- set_names(csv_files, file_path_sans_ext(basename(csv_files))) %>%
  map(read_csv, show_col_types = FALSE)

# -----------------------------
# 1b. Set output directory
# -----------------------------
results_dir <- "outputs/dlmo_resampling_outputs"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 2. Function to resample using PCHIP interpolation
# -----------------------------
resample_profile <- function(df, new_interval_mins) {
  time_numeric <- as.numeric(difftime(df$datetime, min(df$datetime), units = "mins"))
  time_resampled <- seq(0, max(time_numeric), by = new_interval_mins)
  mel_resampled <- pchip(time_numeric, df$melatonin, time_resampled)
  time_resampled_posix <- min(df$datetime) + as.difftime(time_resampled, units = "mins")

  data.frame(datetime = time_resampled_posix, melatonin = mel_resampled)
}

# -----------------------------
# 3. Extract numeric DLMO value
# -----------------------------
extract_dlmo_value <- function(dlmo_result) {
  dlmo_result$ip$inflection_point_fine$x
}

# -----------------------------
# 4. Profile-level function with internal resample-level error capture
# -----------------------------
process_profile <- function(profile_id, df, resample_intervals_mins) {
  original_dlmo_full <- calculate_dlmo(df, threshold = 5)
  original_dlmo_val <- extract_dlmo_value(original_dlmo_full)

  # Safely handle a single resample interval
  safe_inner <- safely(function(interval) {
    resampled_df <- resample_profile(df, interval)
    resampled_dlmo_full <- calculate_dlmo(resampled_df, threshold = 5)
    resampled_dlmo_val <- extract_dlmo_value(resampled_dlmo_full)

    tibble(
      profile = profile_id,
      resample_interval = interval,
      original_dlmo = original_dlmo_val,
      resampled_dlmo = resampled_dlmo_val,
      original_dlmo_full = list(original_dlmo_full),
      resampled_dlmo_full = list(resampled_dlmo_full)
    )
  })

  # Run all intervals in parallel
  results <- future_map(resample_intervals_mins, function(interval) {
    safe_inner(interval)
  })

  # Extract successes
  success <- compact(map(results, "result"))

  # Extract interval-level errors
  errors <- imap_dfr(results, function(res, i) {
    if (!is.null(res$error)) {
      tibble(
        profile = profile_id,
        resample_interval = resample_intervals_mins[[i]],
        error_message = res$error$message
      )
    }
  })

  list(success = bind_rows(success), errors = errors)
}

# -----------------------------
# 5. Setup for parallel run
# -----------------------------
plan(multisession, workers = parallel::detectCores() - 1)
handlers(global = TRUE)

#resample_intervals <- c(2, 5, 10, 15, 20, 30, 45, 60, 75, 90)  # in minutes
resample_intervals <- c(2, 5, 10)  # in minutes

safe_process_profile <- safely(process_profile)

# -----------------------------
# 6. Run DLMO estimation
# -----------------------------
profiles <- head(profiles, 2)
with_progress({
  p <- progressor(along = profiles)
  results_list <- future_imap(profiles, function(df, id) {
    result <- safe_process_profile(id, df, resample_intervals)
    p()
    result
  })
})

# -----------------------------
# 7. Separate and merge results
# -----------------------------
# Successes
all_results <- results_list %>%
  map("result") %>%
  compact() %>%
  map("success") %>%
  bind_rows()

# Interval-level errors
interval_errors <- results_list %>%
  map("result") %>%
  compact() %>%
  map("errors") %>%
  bind_rows()

# Profile-level errors
profile_errors <- imap_dfr(results_list, function(res, id) {
  if (!is.null(res$error)) {
    tibble(profile = id, error_message = res$error$message)
  }
})

if (nrow(all_results) == 0) {
  warning("No successful results to plot. Check errors.")
} else {
  print("all_results contains:")
  print(names(all_results))
  print(head(all_results))
}


# -----------------------------
# 7b. Safely handle empty error tables before merging
# -----------------------------

# Ensure profile_errors has expected structure
if (nrow(profile_errors) == 0) {
  profile_errors <- tibble(profile = character(),
                           error_message = character()) %>%
    mutate(resample_interval = NA_real_)
} else {
  profile_errors <- profile_errors %>%
    mutate(resample_interval = NA_real_)
}

profile_errors <- profile_errors %>%
  select(profile, resample_interval, error_message)

# Ensure interval_errors has expected structure
if (nrow(interval_errors) == 0) {
  interval_errors <- tibble(profile = character(),
                            resample_interval = numeric(),
                            error_message = character())
} else {
  interval_errors <- interval_errors %>%
    select(profile, resample_interval, error_message)
}

# Merge and save
combined_errors <- bind_rows(profile_errors, interval_errors)

# Save outputs
write_csv(combined_errors, file.path(results_dir, "dlmo_resampling_all_errors.csv"))
write_csv(all_results %>% select(-original_dlmo_full, -resampled_dlmo_full),
          file.path(results_dir, "dlmo_resampling_results_summary.csv"))
saveRDS(all_results, file.path(results_dir, "dlmo_resampling_results_full.rds"))

# -----------------------------
# 8. Plot the results
# -----------------------------
# all_results <- readRDS("~/Documents/Projects/DLMO/savedData/dlmo_resamplinganalysis_results_full.rds") # toggle on if loading saved results from directory
if (nrow(all_results) == 0) {
  warning("No successful results to plot.")
} else {
  ggplot(all_results, aes(x = resample_interval, y = resampled_dlmo, color = profile)) +
    geom_line() +
    geom_point() +
    labs(
      title = "DLMO vs. Sampling Interval",
      x = "Sampling Interval (minutes)",
      y = "Estimated DLMO (decimal hours)",
      color = "Profile"
    ) +
    theme_minimal()
}

#plot difference
all_results <- all_results %>%
  mutate(dlmo_difference = resampled_dlmo - original_dlmo)

ggplot(all_results, aes(x = resample_interval, y = dlmo_difference, color = profile)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_line() +
  geom_point() +
  labs(
    title = "Change in DLMO Estimate vs. Sampling Interval",
    x = "Sampling Interval (minutes)",
    y = "Change in DLMO (hours, relative to original)",
    color = "Profile"
  ) +
  theme_minimal()

