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
plan(list(multisession, sequential))
handlers(global = TRUE)

# ──────────────────────────────────
# 2. LOAD AND FILTER PROFILE FILES
# ──────────────────────────────────

profile_folder <- "inst/extdata"


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
#percentages <- c(10, 20)

  scenario2 <- map_dfr(percentages, function(pct) {
    n_del <- floor(pct / 100 * nrow(df))
    if (n_del < 1) return(NULL)

    future_map_dfr(1:50, function(rep) { # discuss reps with Manuel TODO
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
        full_dlmo_time = full_dlmo_time,
        deleted_minutes_from_dlmo = list(
          relative_minutes_to_dlmo(df$datetime[idx], full_dlmo_time)
        )
      )
    }, .options = furrr_options(seed = TRUE))
  })

  return(scenario2)
}

# ────────────────────────────────────────────────────────────────
# 4. RUN ANALYSIS ONLY FOR NEW PROFILES
# ────────────────────────────────────────────────────────────────

results_dir <- "outputs/multiple_deletion_results"
dir.create(results_dir, showWarnings = FALSE)

# Skip already-processed profiles
completed_ids <- list.files(results_dir, pattern = "\\.rds$") %>%
  tools::file_path_sans_ext()
profiles_to_run <- profiles[!names(profiles) %in% completed_ids]

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
num_success <- sum(map_lgl(profiles_to_run, function(id) {
  file.exists(file.path(results_dir, paste0(id, ".rds")))
}))
num_failed <- num_processed - num_success

message("Summary:")
message("Profiles attempted: ", num_processed)
message("  - Successful: ", num_success)
message("  - Failed: ", num_failed)
message("Results saved to '", results_dir, "'")

# Build error log from profiles without .rds files
error_log <- imap_dfr(profiles_to_run, function(df, id) {
  rds_path <- file.path(results_dir, paste0(id, ".rds"))
  if (!file.exists(rds_path)) {
    tibble(profile = id, error_message = "No .rds file generated (likely failed)")
  }
})

# Save error log
write_csv(error_log, file.path(results_dir, "multiple_deletion_errors.csv"))

# ────────────────────────────────────────────────────────────────
# 6. LOAD ALL COMPLETED RESULTS FROM .RDS FILES
# ────────────────────────────────────────────────────────────────

load_all_deletion_results <- function(results_dir) {
  files <- list.files(results_dir, pattern = "\\.rds$", full.names = TRUE)
  valid_files <- files[file.exists(files)]
  safe_read <- purrr::possibly(readRDS, otherwise = NULL)

  purrr::map(valid_files, safe_read) %>%
    purrr::compact() %>%
    dplyr::bind_rows()
}

all_results <- load_all_deletion_results(results_dir)

# Save aggregated results for later reuse
saveRDS(all_results, file.path(results_dir, "dlmo_deletion_all_results.rds"))


# ────────────────────────────────────────────────────────────────
# 7. DEFINE PLOTTING FUNCTIONS
# ────────────────────────────────────────────────────────────────

plot_deletion_raster_by_replicate <- function(results_df) {
  deletions_long <- results_df %>%
    select(profile, percentage_deleted, replicate, deleted_minutes_from_dlmo) %>%
    unnest(deleted_minutes_from_dlmo) %>%
    rename(relative_minute = deleted_minutes_from_dlmo) %>%
    mutate(
      profile = as.factor(profile),
      replicate_id = forcats::fct_inorder(interaction(profile, replicate, sep = ":"))
    )

  ggplot(deletions_long, aes(x = relative_minute, y = replicate_id, fill = profile)) +
    geom_tile(height = 0.8) +
    facet_wrap(~percentage_deleted, scales = "free_y", ncol = 1) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    scale_x_continuous(breaks = seq(-800, 800, by = 60)) +
    labs(
      title = "Deleted Timepoints by Replicate (Relative to DLMO)",
      x = "Time relative to DLMO (minutes)",
      y = "Replicate",
      fill = "Profile"
    ) +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"))
}

plot_dlmo_deletion_violin <- function(results_df) {
  library(gghalves)

  results_df %>%
    filter(scenario == "random_multi") %>%
    ggplot(aes(x = factor(percentage_deleted), y = delta_dlmo)) +
    geom_half_violin(
      aes(fill = factor(percentage_deleted)), side = "l",
      alpha = 0.6, width = 0.9, scale = "width", trim = TRUE,
      bw = 0.2, adjust = 0.5, color = NA
    ) +
    geom_half_point(
      side = "r", shape = 21, size = 1.5, stroke = 0.2,
      color = "black", alpha = 0.6, width = 0.2
    ) +
    geom_half_boxplot(
      side = "r", outlier.shape = NA, width = 0.2,
      color = "black", fill = NA
    ) +
    stat_summary(
      fun = mean, geom = "point", shape = 21,
      size = 2.5, fill = "white", color = "black"
    ) +
    labs(
      title = "DLMO robustness to random deletions",
      subtitle = "Δ DLMO by % of deleted timepoints",
      x = "% Deleted", y = "Δ DLMO (hours)"
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none")
}

# ────────────────────────────────────────────────────────────────
# 8. GENERATE PLOTS
# ────────────────────────────────────────────────────────────────

plot_deletion_raster_by_replicate(all_results)
plot_dlmo_deletion_violin(all_results)
