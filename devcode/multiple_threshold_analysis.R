#-----------------
# load libraries
#-----------------
library(purrr)
library(dplyr)
library(readr)
library(tibble)

# ------------------------------------------------------------------------------
# run_dlmo_across_thresholds()
#
# Wrapper to run DLMO estimation across multiple thresholds for multiple profiles.
#
# For each profile in a named list, this function:
# - Iterates over a user-defined set of thresholds (default: 2, 3, 4, 5, 10)
# - Calls `calculate_dlmo(profile_data, threshold = x)`
# - Extracts:
#     - The full output object (returned by `calculate_dlmo`)
#     - The numeric DLMO estimate (`output$dlmo`, if present)
#     - Any error message (captured via `tryCatch`)
#     - The threshold used
#
# Output files:
# - <results_dir>/full_outputs/<profile>_thresh-<threshold>.rds:
#     - Individual .rds file for each profile × threshold run, storing the full output
# - <results_dir>/dlmo_results_summary.csv:
#     - A clean summary table with one row per run, excluding the heavy `output` column
# - <results_dir>/dlmo_results_combined.rds:
#     - A full tibble of all results including `output`, `dlmo`, `error`, and metadata
#
# Inputs:
# - profiles: named list of melatonin profile data frames
# - thresholds: numeric vector of DLMO thresholds to apply
# - results_dir: directory where output files will be saved
#
# Returns:
# - A tibble with one row per profile × threshold combination, including:
#     - `profile`, `threshold`, `dlmo`, `error`, and the full `output`
# ------------------------------------------------------------------------------
# -------------------------------------------
# Helper to extract numeric DLMO from output
# -------------------------------------------
extract_dlmo_value <- function(dlmo_result) {
  dlmo_result$ip$inflection_point_fine$x
}


run_dlmo_across_thresholds <- function(profiles,
                                       thresholds = c(2, 3, 4, 5, 10),
                                       results_dir = "dlmo_results_thresholds") {
  # Create output folders
  dir.create(file.path(results_dir, "full_outputs"), recursive = TRUE, showWarnings = FALSE)

  # Main loop over profiles
  results <- purrr::map_dfr(names(profiles), function(profile_id) {
    df <- profiles[[profile_id]]

    # Loop over thresholds for each profile
    purrr::map_dfr(thresholds, function(thresh) {
      tryCatch({
        # Run DLMO calculation
        output <- calculate_dlmo(df, threshold = thresh)
        dlmo_val <- tryCatch(
          extract_dlmo_value(output),
          error = function(e) NA_real_
        )

        # Save individual .rds for this run
        saveRDS(
          list(profile = profile_id, threshold = thresh, output = output),
          file = file.path(results_dir, "full_outputs", sprintf("%s_thresh-%s.rds", profile_id, thresh))
        )

        # Return result row
        tibble::tibble(
          profile = profile_id,
          threshold = thresh,
          dlmo = dlmo_val,
          error = NA_character_,
          output = list(output)
        )
      }, error = function(e) {
        # Handle and record any errors
        tibble::tibble(
          profile = profile_id,
          threshold = thresh,
          dlmo = NA_real_,
          error = e$message,
          output = list(NULL)
        )
      })
    })
  })

  # Save clean summary (no list-columns)
  readr::write_csv(
    dplyr::select(results, -output),
    file.path(results_dir, "dlmo_results_summary.csv")
  )

  # Save full results for downstream re-use
  saveRDS(results, file.path(results_dir, "dlmo_results_combined.rds"))

  return(results)
}
#-----------------
# Load profiles
#-----------------
profile_dir <- "inst/extdata"
profile_files <- list.files(profile_dir, pattern = "\\.csv$", full.names = TRUE)
profiles <- profile_files %>%
  set_names(tools::file_path_sans_ext(basename(.))) %>%
  map(read_csv, show_col_types = FALSE)

#------------------------
# Set parameters and run
#------------------------
thresholds_to_test <- c(2, 3, 4, 5, 10)
results_dir <- "outputs/dlmo_threshold_sensitivity"

results <- run_dlmo_across_thresholds(
  profiles = profiles[1:2],
  thresholds = thresholds_to_test,
  results_dir = results_dir
)

#-----------------
# Plot results
#-----------------
library(ggplot2)
library(gghalves)


if (!exists("results")) stop("No 'results' object found. Run analysis first.")
# Filter out failed runs (e.g., those with NA DLMO)
plot_df <- results %>%
  filter(!is.na(dlmo))

# Raincloud-style plot
results_violin<- ggplot(plot_df, aes(x = factor(threshold), y = dlmo)) +
  geom_half_violin(
    aes(fill = factor(threshold)),
    side = "l", alpha = 0.6,
    width = 0.8, scale = "width", trim = TRUE,
    bw = 0.2, adjust = 0.5, color = NA
  ) +
  geom_half_point(
    side = "r", shape = 21, size = 1.5,
    stroke = 0.2, color = "black", alpha = 0.7,
    width = 0.2
  ) +
  geom_half_boxplot(
    side = "r", outlier.shape = NA,
    width = 0.2, color = "black", fill = NA
  ) +
  stat_summary(
    fun = mean, geom = "point", shape = 21,
    size = 2.5, fill = "white", color = "black"
  ) +
  labs(
    title = "DLMO estimates across thresholds",
    x = "Threshold",
    y = "DLMO (decimal hours)",
    fill = "Threshold"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

print(results_violin)

# Summarize number of successful runs per threshold
success_counts <- results %>%
  mutate(success = !is.na(dlmo)) %>%
  count(threshold, success)

# Plot just successful ones
success_plot <- ggplot(success_counts %>% filter(success), aes(x = factor(threshold), y = n)) +
  geom_col(fill = "steelblue") +
  labs(
    title = "Number of successful DLMO estimates per threshold",
    x = "Threshold",
    y = "Number of successful profiles"
  ) +
  theme_minimal(base_size = 14)
print(success_plot)

