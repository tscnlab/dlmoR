# Required packages
library(readr)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pracma)   # for pchip interpolation
library(tools)    # for file name handling
library(furrr)    # for parallelizing
library(future)  # for parallelizing
# -----------------------------
# 1. Load all melatonin profiles from CSVs
# -----------------------------
profile_folder <- system.file("extdata/", package = "dlmoR")  # <- REPLACE THIS with your actual folder path

csv_files <- list.files(profile_folder, pattern = "\\.csv$", full.names = TRUE)

profiles <- set_names(csv_files, file_path_sans_ext(basename(csv_files))) %>%
  map(read_csv, show_col_types = FALSE)

# -----------------------------
# 3. Function to resample using PCHIP interpolation
# -----------------------------
resample_profile <- function(df, new_interval_mins) {
  time_numeric <- as.numeric(difftime(df$datetime, min(df$datetime), units = "mins"))
  time_resampled <- seq(0, max(time_numeric), by = new_interval_mins)
  mel_resampled <- pchip(time_numeric, df$melatonin, time_resampled)
  time_resampled_posix <- min(df$datetime) + as.difftime(time_resampled, units = "mins")

  data.frame(datetime = time_resampled_posix, melatonin = mel_resampled)
}


# -----------------------------
# Function to extract DLMO value from full result
# -----------------------------
extract_dlmo_value <- function(dlmo_result) {
  dlmo_result$ip$inflection_point_fine$x
}

# -----------------------------
# 4. Run DLMO analysis for each profile and resample frequency
# -----------------------------
process_profile <- function(profile_id, df, resample_intervals_mins) {
  original_dlmo_full <- calculate_dlmo(df, threshold = 5)
  original_dlmo_val <- extract_dlmo_value(original_dlmo_full)

  map_dfr(resample_intervals_mins, function(interval) {
    resampled_df <- resample_profile(df, interval)
    resampled_dlmo_full <- calculate_dlmo(resampled_df, threshold = 5)
    resampled_dlmo_val <- extract_dlmo_value(resampled_dlmo_full)

    tibble(
      profile = profile_id,
      resample_interval = interval,
      original_dlmo = original_dlmo_val,
      resampled_dlmo = resampled_dlmo_val,
      original_dlmo_full = list(original_dlmo_full),  # keep full output as list-column
      resampled_dlmo_full = list(resampled_dlmo_full)
    )
  })
}

# -----------------------------
# 5. Run everything
# -----------------------------
# resample_intervals <- c(2, 5, 10, 15, 20, 30, 60, 75)  # in minutes
#
# all_results <- imap_dfr(profiles, function(df, id) {
#   process_profile(id, df, resample_intervals)
# })

# Set parallel backend (multisession works across OSs)
plan(multisession, workers = parallel::detectCores() - 1)  # or specify manually

resample_intervals <- c(2, 5, 10, 15, 20,30 , 45 , 60 , 75 , 90)  # in minutes

# all_results <- imap_dfr(profiles, function(df, id) {
#   process_profile(id, df, resample_intervals)
# })

# Wrap safely
safe_process_profile <- safely(process_profile)

# # Run safely across all profiles # this worked for skipping bad profiles
# results_list <- imap(profiles[13:16], function(df, id) {
#   safe_process_profile(id, df, resample_intervals)
# })

# Parallel version of imap()
results_list <- future_imap(profiles, function(df, id) {
  safe_process_profile(id, df, resample_intervals)
})

# Successful DLMO results
all_results <- compact(map(results_list, "result")) %>% bind_rows()

# Log failed ones
failed_profiles <- imap_dfr(results_list, function(res, id) {
  if (!is.null(res$error)) {
    tibble(profile = id, error_message = res$error$message)
  }
})

write_csv(failed_profiles, "failed_dlmo_profiles.csv")
# -----------------------------
# 6. Plot the results
# -----------------------------
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
