# ------------------------------------------------------------------------------
# DLMO Deletion Robustness Analysis
#
# This script assesses the robustness of DLMO (Dim Light Melatonin Onset)
# estimates to randomly missing data. For each melatonin profile, the script
# computes a baseline DLMO, then deletes 10–50% of timepoints at random (multiple
# replicates per percentage). It re-estimates DLMO after each deletion and records
# the deviation from the original.
#
# Key features:
# - Parallel processing with checkpointed per-profile outputs (.rds files)
# - Nested random deletion with DLMO re-estimation across replicates
# - Error handling and progress tracking
# - Plotting functions for visualizing:
#     1. Deleted timepoints relative to DLMO (raster plot)
#     2. DLMO shifts by % deletion (half-violin plot)
#
# Outputs:
# - Individual .rds files per processed profile in `multiple_deletion_results/`
# - Error log CSV: `multiple_deletion_errors.csv`
# - Aggregated results loaded into `all_results`
# - Visual summary plots (raster + violin)
# ------------------------------------------------------------------------------

# ────────────────────────────────────────────────────────────────
# 1. LOAD LIBRARIES AND SETUP PARALLEL BACKEND
# ────────────────────────────────────────────────────────────────
library(readr)
library(dplyr)
library(purrr)
library(tidyr)
library(future)
library(furrr)
library(progressr)
library(ggplot2)

# Set up nested parallelism (outer = multisession, inner = sequential)
# plan(list(multisession, sequential)) #TODO took out 20250717
#plan(list(multicore, sequential), workers = 16)  # Or whatever limit you want (< 64)
plan(multicore, workers = 64)  # safer than multicore on clusters
handlers(global = TRUE)

# ──────────────────────────────────
# 2. LOAD AND FILTER PROFILE FILES
# ──────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
profile_folder <- if (length(args) > 0) args[[1]] else "inst/extdata"
#profile_folder <- system.file("extdata", package = "dlmoR")

profile_files <- list.files(profile_folder, pattern = "\\.csv$", full.names = TRUE)

profiles <- profile_files %>%
  set_names(tools::file_path_sans_ext(basename(.))) %>%
  map(read_csv, show_col_types = FALSE)

# ────────────────────────────────────────────────────────────────
# 3. DEFINE CORE ANALYSIS FUNCTIONS
# ────────────────────────────────────────────────────────────────

extract_dlmo_value <- function(dlmo_result) {
  dlmo_result$ip$inflection_point_fine$x
}

relative_minutes_to_dlmo <- function(timestamps, dlmo_time) {
  as.numeric(difftime(timestamps, dlmo_time, units = "mins"))
}

run_multiple_deletion <- function(profile_id, df) {
  full_dlmo_result <- tryCatch({
    calculate_dlmo(df, threshold = 5)
  }, error = function(e) {
    message("Full DLMO failed for ", profile_id, ": ", e$message)
    return(NULL)
  })

  if (is.null(full_dlmo_result)) return(tibble())

  full_dlmo <- extract_dlmo_value(full_dlmo_result)
  full_dlmo_time <- dlmoR::decimal_to_posixct(full_dlmo, df$datetime)

  percentages <- c(10, 20, 30, 40, 50)

  scenario2 <- map_dfr(percentages, function(pct) {
    n_del <- floor(pct / 100 * nrow(df))
    if (n_del < 1) return(NULL)

    map_dfr(1:20, function(rep) {
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
        deleted_timepoints = list(df$datetime[idx]),
        delta_dlmo = dlmo_deleted - full_dlmo,
        full_dlmo_dh = full_dlmo,
        deleted_dlmo_dh = dlmo_deleted,
        full_dlmo_time = full_dlmo_time,
        deleted_dlmo_time = if (!is.na(dlmo_deleted))
          dlmoR::decimal_to_posixct(dlmo_deleted, df$datetime)
        else NA,
        deleted_minutes_from_dlmo = list(
          relative_minutes_to_dlmo(df$datetime[idx], full_dlmo_time)
        ),
        replicate_failed = is.na(dlmo_deleted)
      )
    })
  })

  return(scenario2)
}


# ────────────────────────────────────────────────────────────────
# 4. RUN ANALYSIS ONLY FOR NEW PROFILES
# ────────────────────────────────────────────────────────────────

results_dir <- "outputs/multiple_deletion_results"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Skip already-processed profiles
completed_ids <- list.files(results_dir, pattern = "\\.rds$") %>%
  tools::file_path_sans_ext()
profiles_to_run <- profiles[!names(profiles) %in% completed_ids]
#profiles_to_run <- head(profiles_to_run,2)

# Run and save results
safe_run <- safely(run_multiple_deletion)

with_progress({
  p <- progressor(along = profiles_to_run)

  future_imap(profiles_to_run, function(df, id) {
    result <- safe_run(id, df)
    if (!is.null(result$result)) {
      saveRDS(result$result, file = file.path(results_dir, paste0(id, ".rds")))
    }
    p()
    result
  })
})

# ────────────────────────────────────────────────────────────────
# 5. SUMMARY AND ERROR LOGGING
# ────────────────────────────────────────────────────────────────

# Calculate stats
num_processed <- length(profiles_to_run)
num_success <- sum(map_lgl(names(profiles_to_run), function(id) {
  file.exists(file.path(results_dir, paste0(id, ".rds")))
}))
num_failed <- num_processed - num_success

message("Summary:")
message("Profiles attempted: ", num_processed)
message("  - Successful: ", num_success)
message("  - Failed: ", num_failed)
message("Results saved to '", results_dir, "'")

# ---- 1. Profile-level total failures ----
error_log <- tibble(
  profile = names(profiles_to_run)[!map_lgl(names(profiles_to_run), function(id) {
    file.exists(file.path(results_dir, paste0(id, ".rds")))
  })],
  failure_type = "total",
  details = "No results file created"
)

# ---- 2. Load all completed results ----
load_all_deletion_results <- function(results_dir) {
  files <- list.files(results_dir, pattern = "\\.rds$", full.names = TRUE)
  valid_files <- files[file.exists(files)]
  safe_read <- purrr::possibly(readRDS, otherwise = NULL)

  purrr::map(valid_files, safe_read) %>%
    purrr::compact() %>%
    dplyr::bind_rows()
}
all_results <- load_all_deletion_results(results_dir)

# Save aggregated results
saveRDS(all_results, file.path(results_dir, "dlmo_deletion_all_results.rds"))

# ---- 3. Replicate-level partial failures ----
replicate_failures <- all_results %>%
  group_by(profile, percentage_deleted) %>%
  summarise(
    n_reps = n(),
    n_failed = sum(replicate_failed, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_failed > 0) %>%
  mutate(
    failure_type = "partial",
    details = paste0(n_failed, " of ", n_reps,
                     " replicates failed at ", percentage_deleted, "% deleted")
  ) %>%
  select(profile, failure_type, details)

# ---- 4. Combine and save ----
combined_errors <- bind_rows(error_log, replicate_failures) %>%
  arrange(profile, failure_type)

write_csv(combined_errors, file.path(results_dir, "multiple_deletion_error_report.csv"))

