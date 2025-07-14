# ------------------------------------------------------------------------------
# DLMO Single-Point Deletion Sensitivity Analysis
#
# This script quantifies how estimated DLMO (Dim Light Melatonin Onset) values
# respond to the deletion of individual timepoints in melatonin profiles.
# For each profile:
# - The original DLMO is computed
# - Each timepoint is deleted one at a time
# - DLMO is recomputed and compared to the baseline
#
# Key features:
# - Internal parallelization over timepoints via {furrr}
# - Progress bar support via {progressr}
# - Per-profile outputs saved as individual .rds files
# - Final visualization: heatmap of DLMO shifts vs. deletion time
#
# Outputs:
# - Individual .rds result files saved in `single_deletion_results/`
# - Aggregated results loaded and plotted as a DLMO sensitivity heatmap
# ------------------------------------------------------------------------------


# -----------------------------
# Required packages
# -----------------------------
library(readr)
library(dplyr)
library(purrr)
library(tibble)
library(tools)
library(progressr)
library(furrr)
library(future)
library(lubridate)
library(ggplot2)
library(forcats)
library(scico)

# -----------------------------
# Utility functions
# -----------------------------
extract_dlmo_value <- function(dlmo_result) {
  dlmo_result$ip$inflection_point_fine$x
}

relative_minutes_to_dlmo <- function(timestamps, dlmo_time) {
  as.numeric(difftime(timestamps, dlmo_time, units = "mins"))
}

decimal_to_posixct <- function(decimal_hour, reference_times) {
  start_time <- floor_date(min(reference_times), unit = "day")
  start_time + seconds(decimal_hour * 3600)
}

# -----------------------------
# Scenario 1 analysis only (parallelized internally)
# -----------------------------
run_single_deletion <- function(profile_id, df) {
  full_dlmo_result <- tryCatch({
    calculate_dlmo(df, threshold = 5)
  }, error = function(e) {
    message("Full DLMO failed for ", profile_id, ": ", e$message)
    return(NULL)
  })

  if (is.null(full_dlmo_result)) return(tibble())

  full_dlmo <- extract_dlmo_value(full_dlmo_result)
  full_dlmo_time <- decimal_to_posixct(full_dlmo_result$ip$inflection_point_fine$x, full_dlmo_result$prof$datetime)

  results <- future_map_dfr(seq_len(nrow(df)), function(i) {
    df_deleted <- df[-i, ]
    dlmo_deleted <- tryCatch({
      extract_dlmo_value(calculate_dlmo(df_deleted, threshold = 5))
    }, error = function(e) NA_real_)

    tibble(
      profile = profile_id,
      deleted_time = df$datetime[i],
      delta_minutes_from_dlmo = relative_minutes_to_dlmo(df$datetime[i], full_dlmo_time),
      delta_dlmo = dlmo_deleted - full_dlmo
    )
  }, .options = furrr_options(seed = TRUE))

  return(results)
}

# -----------------------------
# Load and process profiles
# -----------------------------
profile_folder <- system.file("extdata/", package = "dlmoR")
csv_files <- list.files(profile_folder, pattern = "\\.csv$", full.names = TRUE)

# Identify processed profiles
result_files <- list.files("single_deletion_results", pattern = "\\.rds$")
processed_ids <- file_path_sans_ext(basename(result_files))

# Load only unprocessed profiles
profiles <- set_names(csv_files, file_path_sans_ext(basename(csv_files)))
profiles <- profiles[!names(profiles) %in% processed_ids]
profiles <- map(profiles, read_csv, show_col_types = FALSE)

# Parallel plan: nested parallelization
plan(nesting(multisession, multisession))
handlers(global = TRUE)

safe_run <- safely(run_single_deletion)

with_progress({
  p <- progressor(along = profiles)
  dir.create("single_deletion_results", showWarnings = FALSE)
  results_list <- future_imap(profiles, function(df, id) {
    result <- safe_run(id, df)
    if (!is.null(result$result)) {
      saveRDS(result$result, file = file.path("single_deletion_results", paste0(id, ".rds")))
    }
    p()
    result
  })
})

# -----------------------------
# Combine and plot results
# -----------------------------
saved_files <- list.files("single_deletion_results", pattern = "\\.rds$", full.names = TRUE)
all_results <- map_dfr(saved_files, readRDS)

# Prepare data for plotting
profile_names <- unique(all_results$profile)
profiles_to_plot <- head(profile_names, -1)

dlmo_range <- range(all_results$delta_dlmo, na.rm = TRUE)
dlmo_range <- c(-max(abs(dlmo_range)), max(abs(dlmo_range)))

filtered_results <- all_results %>%
  filter(profile %in% profiles_to_plot) %>%
  mutate(
    profile = forcats::fct_rev(factor(profile)),
    rounded_minutes = round(delta_minutes_from_dlmo)
  )

plot1 <- ggplot(data = filtered_results, aes(x = rounded_minutes, y = profile, fill = delta_dlmo)) +
  geom_tile(color = NA, height = 0.5) +
  geom_tile(color = NA, height = 0.9, width = 60) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  scale_fill_scico(
    palette = "vik",
    midpoint = 0,
    limits = dlmo_range,
    name = "Δ DLMO (h)"
  ) +
  scale_x_continuous(
    breaks = seq(-800, 800, by = 60),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Sensitivity of DLMO to Single-Point Deletions",
    x = "Deleted Sample Time (minutes relative to DLMO)",
    y = "Profile"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    panel.grid = element_blank()
  )

print(plot1)
