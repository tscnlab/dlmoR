# run_dlmo_threshold_analysis.R

# ------------------------------------------------------------------------------
# DLMO Threshold Sensitivity Analysis Script
# ------------------------------------------------------------------------------
#
# This script runs DLMO (Dim Light Melatonin Onset) calculations across multiple
# threshold levels for a list of melatonin profiles. For each profile × threshold
# combination, it:
# - Runs the DLMO calculation using `calculate_dlmo()`
# - Extracts the DLMO estimate and any error message
# - Saves the full result as an .rds file
# - Collects and saves a summary CSV and full results RDS
# - Optionally, generates diagnostic plots (violin plot, success count)
#
# Output:
# - <results_dir>/full_outputs/<profile>_thresh-<threshold>.rds
# - <results_dir>/dlmo_results_summary.csv
# - <results_dir>/dlmo_results_combined.rds
# - <results_dir>/plots/*.png (optional)
# ------------------------------------------------------------------------------

# Load libraries
library(purrr)
library(dplyr)
library(readr)
library(tibble)
library(ggplot2)
library(gghalves)
library(tools)
library(furrr)
library(future)
library(progressr)

# Set up parallel processing
plan(multisession, workers = parallel::detectCores() - 1)
handlers(global = TRUE)

# Extractor utility
extract_dlmo_value <- function(dlmo_result) {
  dlmo_result$ip$inflection_point_fine$x
}

# Main wrapper function
run_dlmo_across_thresholds <- function(profiles,
                                       thresholds = c(2, 3, 4, 5, 10),
                                       results_dir = "dlmo_results_thresholds") {
  # Create output folder structure
  dir.create(file.path(results_dir, "full_outputs"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(results_dir, "plots"), recursive = TRUE, showWarnings = FALSE)

  # Safe wrapper around DLMO calculation and saving
  safe_run <- safely(function(profile_id, df, threshold) {
    output <- calculate_dlmo(df, threshold = threshold)
    dlmo <- extract_dlmo_value(output)
    saveRDS(list(profile = profile_id, threshold = threshold, output = output),
            file = file.path(results_dir, "full_outputs", sprintf("%s_thresh-%s.rds", profile_id, threshold)))
    tibble(profile = profile_id, threshold = threshold, dlmo = dlmo, error = NA_character_, output = list(output))
  })

  # Create grid of all profile × threshold combinations
  combo_grid <- expand.grid(
    profile_id = names(profiles),
    threshold = thresholds,
    stringsAsFactors = FALSE
  )

  # Parallelized execution with progress bar
  results <- with_progress({
    p <- progressor(along = 1:nrow(combo_grid))

    future_pmap_dfr(combo_grid, function(profile_id, threshold) {
      df <- profiles[[profile_id]]
      result <- safe_run(profile_id, df, threshold)
      p()

      if (!is.null(result$result)) {
        result$result
      } else {
        tibble(profile = profile_id, threshold = threshold, dlmo = NA_real_, error = result$error$message, output = list(NULL))
      }
    }, .options = furrr_options(seed = TRUE))
  })

  # Save summary and full results
  write_csv(select(results, -output), file.path(results_dir, "dlmo_results_summary.csv"))
  saveRDS(results, file.path(results_dir, "dlmo_results_combined.rds"))

  # Optional plotting
  if (nrow(results) > 0) {
    # Violin plot of DLMO values across thresholds
    p_violin <- results %>%
      filter(!is.na(dlmo)) %>%
      ggplot(aes(x = factor(threshold), y = dlmo, fill = factor(threshold))) +
      geom_half_violin(side = "l", alpha = 0.6, width = 0.9, scale = "width", trim = TRUE) +
      geom_half_point(side = "r", shape = 21, size = 1.5, stroke = 0.2, color = "black", alpha = 0.6) +
      geom_half_boxplot(side = "r", outlier.shape = NA, width = 0.2, color = "black", fill = NA) +
      stat_summary(fun = mean, geom = "point", shape = 21, size = 2.5, fill = "white", color = "black") +
      labs(title = "DLMO Estimates by Threshold", x = "Threshold", y = "DLMO (decimal hours)") +
      theme_minimal(base_size = 13) +
      theme(legend.position = "none")

    ggsave(file.path(results_dir, "plots", "dlmo_by_threshold_violin.png"), p_violin, width = 8, height = 6)

    # Success count per threshold
    p_success <- results %>%
      mutate(success = !is.na(dlmo)) %>%
      count(threshold, success) %>%
      filter(success) %>%
      ggplot(aes(x = factor(threshold), y = n)) +
      geom_col(fill = "steelblue") +
      labs(title = "Successful DLMO Estimates by Threshold", x = "Threshold", y = "# Successful Profiles") +
      theme_minimal(base_size = 13)

    ggsave(file.path(results_dir, "plots", "dlmo_success_count_by_threshold.png"), p_success, width = 8, height = 6)
  }

  return(results)
}

# Load profiles
#profile_folder <- system.file("extdata", package = "dlmoR")
profile_folder <- "/home/docker/inst/extdata"

profile_files <- list.files(profile_folder, pattern = "\\.csv$", full.names = TRUE)
profiles <- profile_files %>%
  set_names(tools::file_path_sans_ext(basename(.))) %>%
  map(read_csv, show_col_types = FALSE)

# Parameters
thresholds_to_test <- c(2, 3, 4, 5, 10)
results_dir <- "outputs/dlmo_threshold_sensitivity"

# Run analysis
results <- run_dlmo_across_thresholds(
  profiles = profiles[1:3],
  thresholds = thresholds_to_test,
  results_dir = results_dir
)
