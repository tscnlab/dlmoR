# -----------------------------
# Required packages
# -----------------------------
library(readr)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pracma)
library(tools)
library(furrr)
library(future)
library(progressr)

# -----------------------------
# Utility functions
# -----------------------------
extract_dlmo_value <- function(dlmo_result) {
  dlmo_result$ip$inflection_point_fine$x
}

relative_minutes_to_dlmo <- function(timestamps, dlmo_time) {
  as.numeric(difftime(timestamps, dlmo_time, units = "mins"))
}

# -----------------------------
# Scenario 1 and 2 analysis
# -----------------------------
run_deletion_analysis <- function(profile_id, df) {
  full_dlmo_result <- calculate_dlmo(df, threshold = 5)
  full_dlmo <- extract_dlmo_value(full_dlmo_result)
  full_dlmo_time <- full_dlmo_result$ip$inflection_point_fine$datetime

  # Scenario 1: Single-point deletions
  scenario1 <- map_dfr(seq_len(nrow(df)), function(i) {
    df_deleted <- df[-i, ]
    dlmo_deleted <- tryCatch({
      extract_dlmo_value(calculate_dlmo(df_deleted, threshold = 5))
    }, error = function(e) NA_real_)

    tibble(
      profile = profile_id,
      scenario = "single_point",
      deleted_time = df$datetime[i],
      delta_minutes_from_dlmo = relative_minutes_to_dlmo(df$datetime[i], full_dlmo_time),
      delta_dlmo = dlmo_deleted - full_dlmo
    )
  })

  # Scenario 2: Random deletions by percentage
  percentages <- c(10, 20, 30, 40, 50)
  scenario2 <- map_dfr(percentages, function(pct) {
    n_del <- floor(pct / 100 * nrow(df))
    if (n_del < 1) return(NULL)

    map_dfr(1:10, function(rep) {
      idx <- sample(seq_len(nrow(df)), n_del, replace = FALSE)
      df_deleted <- df[-idx, ]
      dlmo_deleted <- tryCatch({
        extract_dlmo_value(calculate_dlmo(df_deleted, threshold = 5))
      }, error = function(e) NA_real_)

      tibble(
        profile = profile_id,
        scenario = "random_multi",
        percentage_deleted = pct,
        replicate = rep,
        delta_dlmo = dlmo_deleted - full_dlmo
      )
    })
  })

  bind_rows(scenario1, scenario2)
}

# -----------------------------
# Load and process all profiles
# -----------------------------
profile_folder <- system.file("extdata/", package = "dlmoR")
csv_files <- list.files(profile_folder, pattern = "\\.csv$", full.names = TRUE)

profiles <- set_names(csv_files, file_path_sans_ext(basename(csv_files))) %>%
  map(read_csv, show_col_types = FALSE)

profiles_subset <- profiles[1:5]

plan(multisession, workers = parallel::detectCores() - 1)
handlers(global = TRUE)

safe_run <- safely(run_deletion_analysis)

with_progress({
  p <- progressor(along = profiles_subset)
  results_list <- future_imap(profiles_subset, function(df, id) {
    result <- safe_run(id, df)
    p()
    result
  })
})

# Separate successful and failed
all_results <- results_list %>%
  map("result") %>%
  compact() %>%
  bind_rows()

error_log <- imap_dfr(results_list, function(res, id) {
  if (!is.null(res$error)) {
    tibble(profile = id, error_message = res$error$message)
  }
})

write_csv(error_log, "dlmo_deletion_errors.csv")

# -----------------------------
# Plot Scenario 1 heatmap
# -----------------------------
plot1 <- all_results %>%
  filter(scenario == "single_point") %>%
  ggplot(aes(x = delta_minutes_from_dlmo, y = fct_rev(factor(profile)), fill = delta_dlmo)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Δ DLMO (h)") +
  labs(
    title = "Sensitivity of DLMO to Single-Point Deletions",
    x = "Time of Deleted Sample (minutes relative to full-profile DLMO)",
    y = "Profile"
  ) +
  theme_minimal()

# -----------------------------
# Plot Scenario 2 violin plot
# -----------------------------
plot2 <- all_results %>%
  filter(scenario == "random_multi") %>%
  ggplot(aes(x = factor(percentage_deleted), y = delta_dlmo)) +
  geom_violin(aes(fill = factor(percentage_deleted)), alpha = 0.6, width = 0.9) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2, fill = "white") +
  labs(
    title = "DLMO Robustness to Random Deletions (by % removed)",
    x = "Percentage of Deleted Timepoints",
    y = "Delta DLMO (hours)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Display plots
print(plot1)
print(plot2)
