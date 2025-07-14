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
library(ggforce)
library(gghalves)
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

  # # Scenario 1: Single-point deletions
  # scenario1 <- map_dfr(seq_len(nrow(df)), function(i) {
  #   df_deleted <- df[-i, ]
  #   dlmo_deleted <- tryCatch({
  #     extract_dlmo_value(calculate_dlmo(df_deleted, threshold = 5))
  #   }, error = function(e) NA_real_)
  #
  #   tibble(
  #     profile = profile_id,
  #     scenario = "single_point",
  #     deleted_time = df$datetime[i],
  #     delta_minutes_from_dlmo = relative_minutes_to_dlmo(df$datetime[i], full_dlmo_time),
  #     delta_dlmo = dlmo_deleted - full_dlmo
  #   )
  # })

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
# Skip already processed profiles
completed_ids <- list.files("deletion_partial_results", pattern = "\\.rds$") %>%
  file_path_sans_ext()

profiles <- profiles[!names(profiles) %in% completed_ids]
profiles <- profiles[1:2]
profile_folder <- system.file("extdata/", package = "dlmoR")
csv_files <- list.files(profile_folder, pattern = "\\.csv$", full.names = TRUE)

profiles <- set_names(csv_files, file_path_sans_ext(basename(csv_files))) %>%
  map(read_csv, show_col_types = FALSE)

plan(multicore, workers = parallel::detectCores() - 1)
handlers(global = TRUE)

safe_run <- safely(run_deletion_analysis)

with_progress({
  p <- progressor(along = profiles)
  dir.create("deletion_partial_results", showWarnings = FALSE)
  results_list <- future_imap(profiles, function(df, id) {
    result <- safe_run(id, df)
    if (!is.null(result$result)) {
      saveRDS(result$result, file = paste0("deletion_partial_results/", id, ".rds"))
    }
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
# plot1 <- all_results %>%
#   filter(scenario == "single_point") %>%
#   ggplot(aes(x = delta_minutes_from_dlmo, y = factor(profile), fill = delta_dlmo)) +
#   geom_tile(color = "white", linewidth = 0.2) +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Δ DLMO (h)") +
#   labs(
#     title = "Sensitivity of DLMO to Single-Point Deletions",
#     x = "Time of Deleted Sample (minutes relative to full-profile DLMO)",
#     y = "Profile"
#   ) +
#   theme_minimal()

# -----------------------------
# Plot Scenario 2 violin plot
# -----------------------------
# plot2 <- all_results %>%
#   filter(scenario == "random_multi") %>%
#   ggplot(aes(x = factor(percentage_deleted), y = delta_dlmo)) +
#   geom_violin(aes(fill = factor(percentage_deleted)), alpha = 0.6, width = 0.9) +
#   geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
#   stat_summary(fun = mean, geom = "point", shape = 21, size = 2, fill = "white") +
#   labs(
#     title = "DLMO Robustness to Random Deletions (by % removed)",
#     x = "Percentage of Deleted Timepoints",
#     y = "Delta DLMO (hours)"
#   ) +
#   theme_minimal() +
#   theme(legend.position = "none")
#
# # Display plots
# #print(plot1)
# print(plot2)
#
#
# library(ggdist)
#
# plot2 <- all_results %>%
#   filter(scenario == "random_multi") %>%
#   ggplot(aes(x = factor(percentage_deleted), y = delta_dlmo)) +
#
#   # Half-violin (left side)
#   stat_halfeye(
#     aes(fill = factor(percentage_deleted)),
#     adjust = 0.5,
#     width = 0.7,
#     justification = -0.3,  # Push to the left
#     .width = 0,            # Hide interval box
#     point_colour = NA,
#     alpha = 0.6
#   ) +
#
#   # Jittered raw points
#   geom_jitter(width = 0.1, alpha = 0.4, size = 1) +
#
#   # Mean points
#   stat_summary(
#     fun = mean, geom = "point",
#     shape = 21, size = 2, fill = "white"
#   ) +
#
#   labs(
#     title = "DLMO Robustness to Random Deletions (by % removed)",
#     x = "Percentage of Deleted Timepoints",
#     y = "Δ DLMO (hours)"
#   ) +
#   theme_minimal() +
#   theme(legend.position = "none")
#
# print(plot2)
#
# # Install ggforce if needed
# install.packages("ggforce")
# install.packages("gghalves")
# library(gghalves)
#
# plot2 <- all_results %>%
#   filter(scenario == "random_multi") %>%
#   ggplot(aes(x = factor(percentage_deleted), y = delta_dlmo)) +
#
#   # Half violins (left side)
#   geom_half_violin(
#     aes(fill = factor(percentage_deleted)),
#     side = "l",
#     alpha = 0.6,
#     width = 0.8,
#     trim = FALSE
#   ) +
#
#   # Points (right side)
#   geom_jitter(width = 0.15, alpha = 0.5, size = 1, color = "black") +
#
#   # Summary means
#   stat_summary(fun = mean, geom = "point", shape = 21, size = 2, fill = "white") +
#
#   labs(
#     title = "DLMO Robustness to Random Deletions (by % removed)",
#     x = "Percentage of Deleted Timepoints",
#     y = "Δ DLMO (hours)"
#   ) +
#   theme_minimal(base_size = 13) +
#   theme(legend.position = "none")
#
# print(plot2)
#
# library(ggforce)
#
# plot2 <- all_results %>%
#   filter(scenario == "random_multi") %>%
#   ggplot(aes(x = factor(percentage_deleted), y = delta_dlmo)) +
#
#   # Half violins
#   ggforce::geom_sina(aes(color = factor(percentage_deleted)),
#                      alpha = 0.4, size = 0.8, maxwidth = 0.6) +
#
#   ggforce::geom_violinhalf(aes(fill = factor(percentage_deleted)),
#                            side = "l", alpha = 0.6, width = 0.8,
#                            draw_quantiles = c(0.25, 0.5, 0.75)) +
#
#   stat_summary(fun = mean, geom = "point", shape = 21, size = 2, fill = "white") +
#   labs(
#     title = "DLMO Robustness to Random Deletions (by % removed)",
#     x = "Percentage of Deleted Timepoints",
#     y = "Δ DLMO (hours)"
#   ) +
#   theme_minimal(base_size = 13) +
#   theme(legend.position = "none")
#
# print(plot2)
#
# library(ggplot2)
# library(dplyr)
#
# plot2 <- all_results %>%
#   filter(scenario == "random_multi") %>%
#   ggplot(aes(x = factor(percentage_deleted), y = delta_dlmo)) +
#
#   # Full violin plots (mirrored)
#   geom_violin(
#     aes(fill = factor(percentage_deleted)),
#     color = NA,
#     alpha = 0.7,
#     width = 1,
#     scale = "width",
#     trim = TRUE
#   ) +
#
#   # Jittered points (overlaid and centered)
#   geom_jitter(
#     width = 0.15,
#     alpha = 0.5,
#     size = 1.2,
#     shape = 16,
#     color = "black"
#   ) +
#
#   # Mean point (optional)
#   stat_summary(
#     fun = mean,
#     geom = "point",
#     shape = 21,
#     size = 3,
#     fill = "white",
#     color = "black"
#   ) +
#
#   labs(
#     title = "DLMO Robustness to Random Deletions",
#     subtitle = "Δ DLMO by % of deleted timepoints",
#     x = "% Deleted",
#     y = "Δ DLMO (hours)"
#   ) +
#   scale_fill_brewer(palette = "Dark2") +
#   theme_minimal(base_size = 13) +
#   theme(
#     legend.position = "none",
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank()
#   )
#
# print(plot2)
#
# library(gghalves)
#
# plot2 <- all_results %>%
#   filter(scenario == "random_multi") %>%
#   ggplot(aes(x = factor(percentage_deleted), y = delta_dlmo)) +
#
#   geom_half_violin(
#     aes(fill = factor(percentage_deleted)),
#     side = "l",          # 'l' for left, 'r' for right
#     alpha = 0.6,
#     width = 0.9,
#     scale = "width",     # Keep this to fix thinness
#     trim = TRUE,
#     bw = 0.2,
#     adjust = .5,
#     color = NA
#   ) +
#
#   geom_jitter(width = 0.15, alpha = 0.5, size = 1, color = "black") +
#   stat_summary(fun = mean, geom = "point", shape = 21, size = 2, fill = "white") +
#
#   labs(
#     title = "DLMO robustness to random deletions",
#     subtitle = "Δ DLMO by % of deleted timepoints",
#     x = "% Deleted",
#     y = "Δ DLMO (hours)"
#   ) +
#   theme_minimal(base_size = 13) +
#   theme(legend.position = "none")
# print(plot2)
#
# #raincloud
library(gghalves)

plot2 <- all_results %>%
  filter(scenario == "random_multi") %>%
  ggplot(aes(x = factor(percentage_deleted), y = delta_dlmo)) +

  # Half violin on the left
  geom_half_violin(
    aes(fill = factor(percentage_deleted)),
    side = "l",
    alpha = 0.6,
    width = 0.9,
    scale = "width",
    trim = TRUE,
    bw = 0.2,
    adjust = 0.5,
    color = NA
  ) +

  # Points on the right
  geom_half_point(
    side = "r",
    shape = 21,
    size = 1.5,
    stroke = 0.2,
    color = "black",
    alpha = 0.6,
    width = 0.2
  ) +
  geom_half_boxplot(
    side = "r",
    outlier.shape = NA,
    width = 0.2,
    color = "black",
    fill = NA
  )+
  # Optional: Add mean point or boxplot
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 21,
    size = 2.5,
    fill = "white",
    color = "black"
  ) +

  labs(
    title = "DLMO robustness to random deletions",
    subtitle = "Δ DLMO by % of deleted timepoints",
    x = "% Deleted",
    y = "Δ DLMO (hours)"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

print(plot2)

# -----------------------------
# Summary of processing
# -----------------------------
saved_files <- list.files("deletion_partial_results", pattern = "\\.rds$", full.names = TRUE)
all_results <- map_dfr(saved_files, readRDS)

num_processed <- length(results_list)
num_success <- sum(map_lgl(results_list, ~ !is.null(.x$result)))
num_failed <- num_processed - num_success

message("Summary:")
message("Profiles attempted: ", num_processed)
message("  - Successful: ", num_success)
message("  - Failed: ", num_failed)
message("Partial results saved to 'deletion_partial_results/'")
