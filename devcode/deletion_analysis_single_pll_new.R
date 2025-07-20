# ------------------------------------------------------------------------------
# DLMO Single-Point Deletion Sensitivity Analysis
#
# This script quantifies how estimated DLMO (Dim Light Melatonin Onset) values
# respond to the deletion of individual timepoints in melatonin profiles.
#
# For each profile:
# - The original DLMO is computed
# - Each timepoint is deleted one at a time
# - DLMO is recomputed and compared to the baseline
#
# Outputs:
# - Individual .rds result files saved in `single_deletion_results/`
# - Aggregated results visualized as a DLMO sensitivity heatmap
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
# Single deletion analysis (parallelized internally)
# -----------------------------
run_single_deletion <- function(profile_id, df) {
  full_dlmo_result <- tryCatch({
    calculate_dlmo(df, threshold = 2.3)
  }, error = function(e) {
    message("Full DLMO failed for ", profile_id, ": ", e$message)
    return(NULL)
  })

  if (is.null(full_dlmo_result)) return(tibble())

  full_dlmo <- extract_dlmo_value(full_dlmo_result)
  full_dlmo_time <- decimal_to_posixct(full_dlmo_result$ip$inflection_point_fine$x, full_dlmo_result$prof$datetime)

  # Inner parallelism (safe and scoped)
  old_plan <- plan()
  plan(multisession, workers = 2)

  results <- future_map_dfr(seq_len(nrow(df)), function(i) {
    df_deleted <- df[-i, ]
    dlmo_deleted <- tryCatch({
      extract_dlmo_value(calculate_dlmo(df_deleted, threshold = 2.3))
    }, error = function(e) NA_real_)

    tibble(
      profile = profile_id,
      deleted_time = df$datetime[i],
      delta_minutes_from_dlmo = relative_minutes_to_dlmo(df$datetime[i], full_dlmo_time),
      delta_dlmo = dlmo_deleted - full_dlmo
    )
  }, .options = furrr_options(seed = TRUE))

  plan(old_plan)
  return(results)
}


# -----------------------------
# Load and process profiles
# -----------------------------
# Folder with profiles (included with the dlmoR package)
profile_folder <- system.file("extdata/", package = "dlmoR")
csv_files <- list.files(profile_folder, pattern = "\\.csv$", full.names = TRUE)

# Determine already processed profiles
result_files <- list.files("single_deletion_results", pattern = "\\.rds$")
processed_ids <- file_path_sans_ext(basename(result_files))

# Filter to unprocessed profiles
profiles <- set_names(csv_files, file_path_sans_ext(basename(csv_files)))
profiles <- profiles[!names(profiles) %in% processed_ids]

# Optional: only process a subset (e.g. first 2)
# profiles <- head(profiles, 3)

# Read profile data
profiles <- map(profiles, read_csv, show_col_types = FALSE)

# Set outer plan: parallel across profiles
plan(multisession, workers = 6)
handlers(global = TRUE)

safe_run <- safely(run_single_deletion)

# -----------------------------
# Run the sensitivity analysis
# -----------------------------
with_progress({
  p <- progressor(along = profiles)
  dir.create("single_deletion_results", showWarnings = FALSE)

  # results_list <- future_imap(profiles, function(df, id) {
  #   result <- safe_run(id, df)
  #   if (!is.null(result$result)) {
  #     saveRDS(result$result, file = file.path("single_deletion_results", paste0(id, ".rds")))
  #   }
  #   p()
  #   result
  # })
  results_list <- future_imap(profiles, function(df, id) {
    result <- safe_run(id, df)

    if (inherits(result$error, "error")) {
      message(glue::glue("Profile {id} failed: {result$error$message}"))
    } else {
      message(glue::glue("Profile {id} succeeded."))
      saveRDS(result$result, file = file.path("single_deletion_results", paste0(id, ".rds")))
    }

    p()
    result
  })

})

# -----------------------------
# Combine and plot results
# -----------------------------
# Read individual result files
saved_files <- list.files("single_deletion_results", pattern = "\\.rds$", full.names = TRUE)
all_results <- map_dfr(saved_files, readRDS)

# Prepare data for plotting
profile_names <- unique(all_results$profile)
profiles_to_plot <- head(profile_names, -1)

dlmo_range <- range(all_results$delta_dlmo, na.rm = TRUE)
dlmo_range <- c(-max(abs(dlmo_range)), max(abs(dlmo_range)))  # symmetric limits
#
# filtered_results <- all_results %>%
#   filter(profile %in% profiles_to_plot) %>%
#   mutate(
#     profile = forcats::fct_rev(factor(profile)),
#     rounded_minutes = round(delta_minutes_from_dlmo)
#   )

# # Plot DLMO sensitivity heatmap
# plot1 <- ggplot(filtered_results, aes(x = rounded_minutes, y = profile, fill = delta_dlmo)) +
#   geom_tile(color = NA, height = 0.5) +
#   geom_tile(color = NA, height = 0.9, width = 60) +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
#   scale_fill_scico(
#     palette = "vik",
#     midpoint = 0,
#     limits = dlmo_range,
#     name = "Δ DLMO (h)"
#   ) +
#   scale_x_continuous(
#     breaks = seq(-800, 800, by = 60),
#     expand = expansion(mult = c(0, 0))
#   ) +
#   coord_cartesian(clip = "off") +
#   labs(
#     title = "Sensitivity of DLMO to Single-Point Deletions",
#     x = "Deleted Sample Time (minutes relative to DLMO)",
#     y = "Profile"
#   ) +
#   theme_minimal(base_size = 12) +
#   theme(
#     axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
#     panel.grid = element_blank()
#   )
#
# print(plot1)
#
#
# ##
# # Load required packages
# library(ggplot2)
# library(dplyr)
# library(forcats)
# library(scico)
# library(patchwork)  # for stacking plots
#
# # Assuming `all_results` already loaded from deletion_scenario1_results/*.rds
# # Filter and prepare
# dlmo_range <- range(all_results$delta_dlmo, na.rm = TRUE)
# dlmo_range <- c(-max(abs(dlmo_range)), max(abs(dlmo_range)))
#
# filtered_results <- all_results %>%
#   filter(!is.na(delta_dlmo)) %>%
#   mutate(
#     profile = forcats::fct_rev(factor(profile)),
#     rounded_minutes = round(delta_minutes_from_dlmo)
#   )
#
# # --- Ribbon Plot (Top) ---
# summary_data <- filtered_results %>%
#   mutate(binned_minutes = floor(rounded_minutes / 10) * 10) %>%
#   group_by(binned_minutes) %>%
#   summarise(
#     mean_dlmo = mean(delta_dlmo, na.rm = TRUE),
#     sd_dlmo = sd(delta_dlmo, na.rm = TRUE),
#     .groups = "drop"
#   )
#
# plot_top <- ggplot(summary_data, aes(x = binned_minutes, y = mean_dlmo)) +
#   geom_ribbon(aes(ymin = mean_dlmo - sd_dlmo, ymax = mean_dlmo + sd_dlmo), fill = "skyblue", alpha = 0.4) +
#   geom_line(color = "black", linewidth = 1) +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
#   geom_hline(yintercept = 0, linetype = "dotted", color = "gray60") +
#   labs(
#     y = "Δ DLMO (h)",
#     x = NULL,
#     title = NULL
#   ) +
#   scale_x_continuous(
#     breaks = seq(-800, 800, by = 120),
#     expand = expansion(mult = c(0, 0))
#   ) +
#   theme_minimal(base_size = 13) +
#   theme(
#     axis.text.x = element_blank(),
#     axis.title.x = element_blank(),
#     panel.grid.minor = element_blank()
#   )
#
# # --- Heatmap (Bottom) ---
# plot_bottom <- ggplot(filtered_results, aes(x = rounded_minutes, y = profile, fill = delta_dlmo)) +
#   geom_tile(color = NA, height = 0.5) +
#   geom_tile(color = NA, height = 0.9, width = 60) +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
#   scale_fill_scico(
#     palette = "vik",
#     midpoint = 0,
#     limits = dlmo_range,
#     name = "Δ DLMO (h)"
#   ) +
#   scale_x_continuous(
#     breaks = seq(-800, 800, by = 120),
#     expand = expansion(mult = c(0, 0))
#   ) +
#   coord_cartesian(clip = "off") +
#   labs(
#     x = "Deleted timepoint (minutes relative to DLMO)",
#     y = "Melatonin profile"
#   ) +
#   theme_minimal(base_size = 13) +
#   theme(
#     axis.text.x = element_text(angle = 0, vjust = 0.5),
#     axis.text.y = element_blank(),
#     axis.ticks.y = element_blank(),
#     panel.grid = element_blank()
#   )
#
# # --- Combine with patchwork ---
# stacked_plot <- plot_top / plot_bottom + plot_layout(heights = c(1, 2))
# print(stacked_plot)


library(dplyr)
library(ggplot2)
library(forcats)
library(scico)
library(patchwork)

# Use binned minutes for alignment
filtered_results <- all_results %>%
  filter(!is.na(delta_dlmo)) %>%
  mutate(
    profile = forcats::fct_rev(factor(profile)),
    binned_minutes = floor(delta_minutes_from_dlmo / 10) * 10  # consistent binning
  )

# Top: mean + SD ribbon summary
summary_data <- filtered_results %>%
  group_by(binned_minutes) %>%
  summarise(
    mean_dlmo = mean(delta_dlmo, na.rm = TRUE),
    sd_dlmo = sd(delta_dlmo, na.rm = TRUE),
    n = sum(!is.na(delta_dlmo)),
    sem_dlmo = sd_dlmo / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(group = "Bias (Mean ± SD)")

# Shared x-axis scale
x_breaks <- seq(min(filtered_results$binned_minutes), max(filtered_results$binned_minutes), by = 60)
shared_x <- scale_x_continuous(
  breaks = x_breaks,
  expand = expansion(mult = c(0, 0))
)

# Top plot: Ribbon and smoothed line
plot_top <- ggplot(summary_data, aes(x = binned_minutes)) +
  geom_ribbon(
    aes(ymin = mean_dlmo - sd_dlmo, ymax = mean_dlmo + sd_dlmo, fill = group),
    alpha = 0.4
  ) +
  geom_line(
    aes(y = mean_dlmo, color = group),
    linewidth = 1
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray60") +
  shared_x +
  scale_fill_manual(name = NULL, values = c("Bias (Mean ± SD)" = "skyblue")) +
  scale_color_manual(name = NULL, values = c("Bias (Mean ± SD)" = "black")) +
  labs(
    y = "Δ DLMO (h)",
    x = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top",
    legend.justification = "left"
  )

# Bottom plot: Heatmap of DLMO shifts
dlmo_range <- range(filtered_results$delta_dlmo, na.rm = TRUE)
dlmo_range <- c(-max(abs(dlmo_range)), max(abs(dlmo_range)))  # symmetric

plot_bottom <- ggplot(filtered_results, aes(x = binned_minutes, y = profile, fill = delta_dlmo)) +
  geom_tile(color = NA, height = 0.9, width = 10) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  shared_x +
  scale_fill_scico(
    palette = "vik",
    midpoint = 0,
    limits = dlmo_range,
    name = "Δ DLMO (h)"
  ) +
  labs(
    x = "Deleted timepoint (minutes relative to DLMO)",
    y = "Melatonin profile"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank()
  )

# Combine
plot_stacked <- plot_top / plot_bottom + plot_layout(heights = c(1, 2))
print(plot_stacked)
