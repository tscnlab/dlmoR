# Required packages
library(readr)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pracma)   # for pchip interpolation
library(tools)    # for file name handling
library(furrr)    # for parallelizing
library(future)   # for parallelizing

#-----------------------------
# Helper Functions
#-----------------------------
resample_profile <- function(df, new_interval_mins) {
  time_numeric <- as.numeric(difftime(df$datetime, min(df$datetime), units = "mins"))
  time_resampled <- seq(0, max(time_numeric), by = new_interval_mins)
  mel_resampled <- pchip(time_numeric, df$melatonin, time_resampled)
  time_resampled_posix <- min(df$datetime) + as.difftime(time_resampled, units = "mins")
  data.frame(datetime = time_resampled_posix, melatonin = mel_resampled)
}

extract_dlmo_value <- function(dlmo_result) {
  dlmo_result$ip$inflection_point_fine$x
}

#-----------------------------
# Main Deletion Analysis Function
#-----------------------------
run_deletion_analysis <- function(profile_id, df) {
  full_dlmo_result <- calculate_dlmo(df, threshold = 5)
  full_dlmo <- extract_dlmo_value(full_dlmo_result)

  # Scenario 1: Single-point deletion
  scenario1_results <- map_dfr(seq_len(nrow(df)), function(i) {
    df_deleted <- df[-i, ]
    dlmo_deleted <- tryCatch({
      extract_dlmo_value(calculate_dlmo(df_deleted, threshold = 5))
    }, error = function(e) NA_real_)

    tibble(
      profile = profile_id,
      scenario = "single_point",
      deleted_index = i,
      deleted_time = df$datetime[i],
      delta_dlmo = dlmo_deleted - full_dlmo
    )
  })

  # Scenario 2: Random deletions with varying N
  max_N <- floor(nrow(df) * 0.5)
  scenario2_results <- map_dfr(2:max_N, function(n_del) {
    map_dfr(1:10, function(rep) {
      deleted_idx <- sample(seq_len(nrow(df)), n_del, replace = FALSE)
      df_deleted <- df[-deleted_idx, ]
      dlmo_deleted <- tryCatch({
        extract_dlmo_value(calculate_dlmo(df_deleted, threshold = 5))
      }, error = function(e) NA_real_)

      tibble(
        profile = profile_id,
        scenario = "random_multi",
        n_deleted = n_del,
        replicate = rep,
        delta_dlmo = dlmo_deleted - full_dlmo
      )
    })
  })

  bind_rows(scenario1_results, scenario2_results)
}

#-----------------------------
# Run All Profiles in Parallel
#-----------------------------
# Replace with actual profile folder
profile_folder <- system.file("extdata/", package = "dlmoR")
csv_files <- list.files(profile_folder, pattern = "\\.csv$", full.names = TRUE)

profiles <- set_names(csv_files, file_path_sans_ext(basename(csv_files))) %>%
  map(read_csv, show_col_types = FALSE)

plan(multisession, workers = parallel::detectCores() - 1)

results_list <- future_imap(profiles, function(df, id) {
  safely(run_deletion_analysis)(id, df)
})

all_results <- compact(map(results_list, "result")) %>% bind_rows()

#-----------------------------
# Plot Scenario 1: Single-point Deletion
#-----------------------------
plot1 <- all_results %>%
  filter(scenario == "single_point") %>%
  ggplot(aes(x = deleted_time, y = delta_dlmo, color = profile)) +
  geom_point() +
  geom_line() +
  labs(
    title = "Change in DLMO from Single-Point Deletions",
    x = "Time of Deleted Sample",
    y = "Delta DLMO (hours)"
  ) +
  theme_minimal()

#-----------------------------
# Plot Scenario 2: Random Multi-Point Deletions
#-----------------------------
plot2 <- all_results %>%
  filter(scenario == "random_multi") %>%
  ggplot(aes(x = factor(n_deleted), y = delta_dlmo)) +
  geom_violin(aes(fill = factor(n_deleted)), alpha = 0.6, width = 0.9) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2, fill = "white") +
  labs(
    title = "DLMO Robustness to Random Multi-Point Deletions",
    x = "Number of Deleted Timepoints (N)",
    y = "Delta DLMO (hours)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Display plots
print(plot1)
print(plot2)
