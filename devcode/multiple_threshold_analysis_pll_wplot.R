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

  # Run analysis with safe error capture
  safe_run <- safely(function(profile_id, df, threshold) {
    output <- calculate_dlmo(df, threshold = threshold)
    dlmo <- extract_dlmo_value(output)
    saveRDS(list(profile = profile_id, threshold = threshold, output = output),
            file = file.path(results_dir, "full_outputs", sprintf("%s_thresh-%s.rds", profile_id, threshold)))
    tibble(profile = profile_id, threshold = threshold, dlmo = dlmo, error = NA_character_, output = list(output))
  })

  # Parallelized with progress bar
  results <- with_progress({
    p <- progressor(along = names(profiles))

    future_map_dfr(names(profiles), function(profile_id) {
      df <- profiles[[profile_id]]

      map_dfr(thresholds, function(thresh) {
        result <- safe_run(profile_id, df, thresh)

        if (!is.null(result$result)) {
          result$result
        } else {
          tibble(profile = profile_id, threshold = thresh, dlmo = NA_real_, error = result$error$message, output = list(NULL))
        }
      }) %>% { p(); . }

    }, .options = furrr_options(seed = TRUE))
  })

  # Save summary and full results
  write_csv(select(results, -output), file.path(results_dir, "dlmo_results_summary.csv"))
  saveRDS(results, file.path(results_dir, "dlmo_results_combined.rds"))

  # # Optional plotting
  # if (nrow(results) > 0) {
  #   # Violin plot of DLMO values across thresholds
  #   p_violin <- results %>%
  #     filter(!is.na(dlmo)) %>%
  #     ggplot(aes(x = factor(threshold), y = dlmo, fill = factor(threshold))) +
  #     geom_half_violin(side = "l", alpha = 0.6, width = 0.9, scale = "width", trim = TRUE) +
  #     geom_half_point(side = "r", shape = 21, size = 1.5, stroke = 0.2, color = "black", alpha = 0.6) +
  #     geom_half_boxplot(side = "r", outlier.shape = NA, width = 0.2, color = "black", fill = NA) +
  #     stat_summary(fun = mean, geom = "point", shape = 21, size = 2.5, fill = "white", color = "black") +
  #     labs(title = "DLMO Estimates by Threshold", x = "Threshold", y = "DLMO (decimal hours)") +
  #     theme_minimal(base_size = 13) +
  #     theme(legend.position = "none")
  #
  #   ggsave(file.path(results_dir, "plots", "dlmo_by_threshold_violin.png"), p_violin, width = 8, height = 6)
  #
  #   # Success count per threshold
  #   p_success <- results %>%
  #     mutate(success = !is.na(dlmo)) %>%
  #     count(threshold, success) %>%
  #     filter(success) %>%
  #     ggplot(aes(x = factor(threshold), y = n)) +
  #     geom_col(fill = "steelblue") +
  #     labs(title = "Successful DLMO Estimates by Threshold", x = "Threshold", y = "# Successful Profiles") +
  #     theme_minimal(base_size = 13)
  #
  #   ggsave(file.path(results_dir, "plots", "dlmo_success_count_by_threshold.png"), p_success, width = 8, height = 6)
  # }

  return(results)
}

# Load profiles
profile_dir <- "inst/extdata"
profile_files <- list.files(profile_dir, pattern = "\\.csv$", full.names = TRUE)
profiles <- profile_files %>%
  set_names(tools::file_path_sans_ext(basename(.))) %>%
  map(read_csv, show_col_types = FALSE)

# Parameters
thresholds_to_test <- c(2, 3, 4, 5, 10)
results_dir <- "outputs/dlmo_threshold_sensitivity"

# Run analysis
results <- run_dlmo_across_thresholds(
  profiles = profiles,
  thresholds = thresholds_to_test,
  results_dir = results_dir
)


####################
# post hoc plotting
####################
results <- readRDS("~/Documents/Projects/DLMO/dlmoRpaperresults/Blume/threshold_results/dlmo_threshold_sensitivity/dlmo_results_combined.rds") # toggle on if loading saved results from directory

library(gghalves)
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

  #ggsave(file.path(results_dir, "plots", "dlmo_by_threshold_violin.png"), p_violin, width = 8, height = 6)

  # Success count per threshold
  p_success <- results %>%
    mutate(success = !is.na(dlmo)) %>%
    count(threshold, success) %>%
    filter(success) %>%
    ggplot(aes(x = factor(threshold), y = n)) +
    geom_col(fill = "steelblue") +
    labs(title = "Successful DLMO Estimates by Threshold", x = "Threshold", y = "# Successful Profiles") +
    theme_minimal(base_size = 13)

  #ggsave(file.path(results_dir, "plots", "dlmo_success_count_by_threshold.png"), p_success, width = 8, height = 6)
}


##
# Handle cases where threshold 5 failed by filtering profiles that have threshold 5
results_delta <- results %>%
  filter(!is.na(dlmo)) %>%
  group_by(profile) %>%
  filter(any(threshold == 10 & !is.na(dlmo))) %>%  # keep only profiles where threshold 5 succeeded
  mutate(dlmo_ref = dlmo[threshold == 10][1]) %>%  # take first available ref if multiple
  ungroup() %>%
  mutate(delta_dlmo = dlmo - dlmo_ref)

# Plot change in DLMO vs threshold
library(gghalves)

delta_dlmo_plot<- ggplot(results_delta, aes(x = factor(threshold), y = delta_dlmo, fill = factor(threshold))) +
  geom_half_violin(side = "l", alpha = 0.6, width = 0.9, scale = "width", trim = TRUE) +
  geom_half_point(side = "r", shape = 21, size = 1.5, stroke = 0.2, color = "black", alpha = 0.6) +
  geom_half_boxplot(side = "r", outlier.shape = NA, width = 0.2, color = "black", fill = NA, position = position_nudge(x = 0.4)) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2.5, fill = "white", color = "black", position = position_nudge(x = 0.4)) +
  labs(
    title = "Change in DLMO vs Threshold",
    subtitle = "Relative to DLMO at threshold = 10",
    x = "Threshold",
    y = "Δ DLMO (decimal hours)"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

print(delta_dlmo_plot)

# Count number of successful DLMO estimates at each threshold
library(dplyr)

# summary
results %>%
  filter(!is.na(dlmo)) %>%
  count(threshold, name = "n_successful") %>%
  arrange(threshold)

# total # of profiles
results %>% distinct(profile) %>% count()

# plot % success at each threshld
library(dplyr)
library(ggplot2)

# Compute % success per threshold
success_summary <- results %>%
  group_by(threshold) %>%
  summarise(
    n_total = n(),
    n_success = sum(!is.na(dlmo)),
    pct_success = 100 * n_success / n_total,
    .groups = "drop"
  )
print(success_summary)

# Line plot of % success
psuccess <- ggplot(success_summary, aes(x = threshold, y = pct_success)) +
  geom_line(color = "steelblue", size = 1) +
  geom_point(color = "steelblue", size = 2) +
  scale_x_continuous(breaks = success_summary$threshold) +  # only show ticks at data points
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  labs(
    title = "Success Rate of DLMO Estimates by Threshold",
    x = "Threshold [pg/mL]",
    y = "% of Profiles with Valid DLMO"
  ) +
  theme_minimal(base_size = 14)
print(psuccess)

library(ggplot2)
library(patchwork)


##############
# combo plots
##############
#results <- readRDS("~/Documents/Projects/DLMO/dlmoRpaperresults/Blume/threshold_results/dlmo_threshold_sensitivity/dlmo_results_combined.rds") # toggle on if loading saved results from directory
results <- readRDS("~/Documents/Projects/DLMO/dlmoRpaperresults/Civibe/threshold/outputs/dlmo_threshold_sensitivity/dlmo_results_combined.rds") # toggle on if loading saved results from directory

results_delta <- results %>%
  filter(!is.na(dlmo)) %>%
  group_by(profile) %>%
  filter(any(threshold == 2 & !is.na(dlmo))) %>%  # keep only profiles where threshold 5 succeeded
  mutate(dlmo_ref = dlmo[threshold == 2][1]) %>%  # take first available ref if multiple
  ungroup() %>%
  mutate(delta_dlmo = dlmo - dlmo_ref)

results_delta <- results %>%
  filter(!is.na(dlmo)) %>%
  group_by(profile) %>%
  filter(any(threshold == 2 & !is.na(dlmo))) %>%
  mutate(dlmo_ref = dlmo[threshold == 2][1]) %>%
  ungroup() %>%
  mutate(delta_dlmo = dlmo - dlmo_ref) %>%
  filter(abs(delta_dlmo) <= 4)  # Keep only differences within ±6 hours

# results_delta <- results_delta %>%
#   mutate(delta_dlmo = if_else(threshold == 2, delta_dlmo + 0.0000000000000001, delta_dlmo))

# Your violin plot (e.g., p_violin)
point_nudge <- 0.08

p_delta_dlmo <- ggplot(results_delta, aes(x = factor(threshold), y = delta_dlmo, fill = factor(threshold))) +
  geom_half_violin(side = "l", alpha = 0.6, width = 0.9, scale = "width", trim = TRUE) +
  geom_point(
  aes(x = as.numeric(factor(threshold)) + point_nudge, y = delta_dlmo),
  shape = 21,
  size = 1.5,
  stroke = 0.2,
  color = "black",
  alpha = 0.6,
  position = position_jitter(width = 0.05, height = 0)  # mild horizontal jitter
) +
  geom_half_boxplot(side = "r", outlier.shape = NA, width = 0.2, color = "black", fill = NA, position = position_nudge(x = 0.18)) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2.5, fill = "white", color = "black", position = position_nudge(x = 0.18)) +
  labs(
    title = "Influence of threshold parameter on DLMO estimates",
    subtitle = "Relative to DLMO at threshold = 2 pg/mL",
    x = "Threshold",
    y = "Δ DLMO (decimal hours)"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")


delta_summary <- results_delta %>%
  group_by(threshold) %>%
  summarise(
    mean_delta = mean(delta_dlmo, na.rm = TRUE),
    sd_delta = sd(delta_dlmo, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

print(delta_summary)
write.csv(delta_summary, "civibe_delta_dlmo_summary_threshold.csv", row.names = FALSE)

# Your % success plot (e.g., p_success)
success_summary <- results %>%
  group_by(threshold) %>%
  summarise(n_success = sum(!is.na(dlmo)),
            total = n_distinct(profile),
            pct_success = 100 * n_success / total,
            .groups = "drop")

write.csv(success_summary, "civibe_success_summary_threshold.csv", row.names = FALSE)

p_success <- ggplot(success_summary, aes(x = threshold, y = pct_success)) +
  geom_line(color = "steelblue", size = 1) +
  geom_point(color = "steelblue", size = 2) +
  scale_x_continuous(breaks = success_summary$threshold) +  # only show ticks at data points
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  labs(
    title = "Success rate of DLMO detection",
    x = "Threshold [pg/mL]",
    y = "% of Profiles with Valid DLMO"
  ) +
  theme_minimal(base_size = 14)

# Combine using patchwork with relative widths
combined_plot <- p_delta_dlmo + p_success + plot_layout(widths = c(3, 1))

# Display
print(combined_plot)
