# ------------------------------------------------------------------------------
# DLMO Resampling Sensitivity Analysis
#
# This script evaluates how estimated DLMO (Dim Light Melatonin Onset) values
# change as a function of temporal sampling interval. It loads multiple melatonin
# profiles, computes the original DLMO for each, then resamples each profile at
# various intervals using PCHIP interpolation. DLMO is re-estimated for each
# resampled profile, and the differences are recorded and visualized.
#
# Key features:
# - Parallel processing across profiles and resample intervals (via {furrr})
# - Error handling at both profile and interval levels
# - Progress bar support via {progressr}
# - Output includes results table, error logs, and summary plots
#
# Outputs:
# - RDS file of all successful results
# - CSV file of all errors
# - Two plots:
#     1. DLMO vs. Sampling Interval
#     2. Change in DLMO vs. Sampling Interval
# ------------------------------------------------------------------------------

# -----------------------------
# Required packages
# -----------------------------
library(readr)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pracma)    # for pchip interpolation
library(tools)     # for file name handling
library(furrr)     # for parallelization
library(future)    # for parallelization
library(progressr) # for progress bar

# -----------------------------
# 1a. Load all melatonin profiles from CSVs
# -----------------------------
#profile_folder <- system.file("extdata/", package = "dlmoR")
args <- commandArgs(trailingOnly = TRUE)
profile_folder <- if (length(args) > 0) args[[1]] else "inst/extdata"

csv_files <- list.files(profile_folder, pattern = "\\.csv$", full.names = TRUE)

profiles <- set_names(csv_files, file_path_sans_ext(basename(csv_files))) %>%
  map(read_csv, show_col_types = FALSE)

# -----------------------------
# 1b. Set output directory
# -----------------------------
results_dir <- "outputs/dlmo_resampling_outputs"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 2. Function to resample using PCHIP interpolation
# -----------------------------
resample_profile <- function(df, new_interval_mins) {
  time_numeric <- as.numeric(difftime(df$datetime, min(df$datetime), units = "mins"))
  time_resampled <- seq(0, max(time_numeric), by = new_interval_mins)
  mel_resampled <- pchip(time_numeric, df$melatonin, time_resampled)
  time_resampled_posix <- min(df$datetime) + as.difftime(time_resampled, units = "mins")

  data.frame(datetime = time_resampled_posix, melatonin = mel_resampled)
}

# -----------------------------
# 3. Extract numeric DLMO value
# -----------------------------
extract_dlmo_value <- function(dlmo_result) {
  dlmo_result$ip$inflection_point_fine$x
}

# -----------------------------
# 4. Profile-level function with internal resample-level error capture
# -----------------------------
process_profile <- function(profile_id, df, resample_intervals_mins) {
  original_dlmo_full <- calculate_dlmo(df, threshold = 5)
  original_dlmo_val <- extract_dlmo_value(original_dlmo_full)

  # Safely handle a single resample interval
  safe_inner <- safely(function(interval) {
    resampled_df <- resample_profile(df, interval)
    resampled_dlmo_full <- calculate_dlmo(resampled_df, threshold = 5)
    resampled_dlmo_val <- extract_dlmo_value(resampled_dlmo_full)

    tibble(
      profile = profile_id,
      resample_interval = interval,
      original_dlmo = original_dlmo_val,
      resampled_dlmo = resampled_dlmo_val,
      original_dlmo_full = list(original_dlmo_full),
      resampled_dlmo_full = list(resampled_dlmo_full)
    )
  })

  # Run all intervals in parallel
  results <- future_map(resample_intervals_mins, function(interval) {
    safe_inner(interval)
  })

  # Extract successes
  success <- compact(map(results, "result"))

  # Extract interval-level errors
  errors <- imap_dfr(results, function(res, i) {
    if (!is.null(res$error)) {
      tibble(
        profile = profile_id,
        resample_interval = resample_intervals_mins[[i]],
        error_message = res$error$message
      )
    }
  })

  list(success = bind_rows(success), errors = errors)
}

# -----------------------------
# 5. Setup for parallel run
# -----------------------------
plan(multisession, workers = 16)
handlers(global = TRUE)

resample_intervals <- c(2, 5, 10, 15, 20, 30, 45, 60, 75, 90)  # in minutes
#resample_intervals <- c(2, 5, 10)  # in minutes

safe_process_profile <- safely(process_profile)

# -----------------------------
# 6. Run DLMO estimation
# -----------------------------
#profiles <- head(profiles, 2)
with_progress({
  p <- progressor(along = profiles)
  results_list <- future_imap(profiles, function(df, id) {
    result <- safe_process_profile(id, df, resample_intervals)
    p()
    result
  })
})

# -----------------------------
# 7. Separate and merge results
# -----------------------------
# Successes
all_results <- results_list %>%
  map("result") %>%
  compact() %>%
  map("success") %>%
  bind_rows()

# Interval-level errors
interval_errors <- results_list %>%
  map("result") %>%
  compact() %>%
  map("errors") %>%
  bind_rows()

# Profile-level errors
profile_errors <- imap_dfr(results_list, function(res, id) {
  if (!is.null(res$error)) {
    tibble(profile = id, error_message = res$error$message)
  }
})

if (nrow(all_results) == 0) {
  warning("No successful results to plot. Check errors.")
} else {
  print("all_results contains:")
  print(names(all_results))
  print(head(all_results))
}


# -----------------------------
# 7b. Safely handle empty error tables before merging
# -----------------------------

# Ensure profile_errors has expected structure
if (nrow(profile_errors) == 0) {
  profile_errors <- tibble(profile = character(),
                           error_message = character()) %>%
    mutate(resample_interval = NA_real_)
} else {
  profile_errors <- profile_errors %>%
    mutate(resample_interval = NA_real_)
}

profile_errors <- profile_errors %>%
  select(profile, resample_interval, error_message)

# Ensure interval_errors has expected structure
if (nrow(interval_errors) == 0) {
  interval_errors <- tibble(profile = character(),
                            resample_interval = numeric(),
                            error_message = character())
} else {
  interval_errors <- interval_errors %>%
    select(profile, resample_interval, error_message)
}

# Merge and save
combined_errors <- bind_rows(profile_errors, interval_errors)

# Save outputs
write_csv(combined_errors, file.path(results_dir, "dlmo_resampling_all_errors.csv"))
write_csv(all_results %>% select(-original_dlmo_full, -resampled_dlmo_full),
          file.path(results_dir, "dlmo_resampling_results_summary.csv"))
saveRDS(all_results, file.path(results_dir, "dlmo_resampling_results_full.rds"))

# -----------------------------
# 8. Plot the results
# -----------------------------
# all_results <- readRDS("~/Documents/Projects/DLMO/savedData/dlmo_resamplinganalysis_results_full.rds") # toggle on if loading saved results from directory
if (nrow(all_results) == 0) {
  warning("No successful results to plot.")
} else {
  ggplot(all_results, aes(x = resample_interval, y = resampled_dlmo, color = profile)) +
    geom_line() +
    geom_point() +
    labs(
      title = "DLMO vs. Sampling Interval",
      x = "Sampling Interval (minutes)",
      y = "Estimated DLMO (decimal hours)",
      color = "Profile"
    )+
    theme_minimal()+
    theme(legend.position = "none")
}

#plot difference
all_results <- all_results %>%
  mutate(dlmo_difference = resampled_dlmo - original_dlmo)

ggplot(all_results, aes(x = resample_interval, y = dlmo_difference, color = profile)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_line() +
  geom_point() +
  labs(
    title = "Change in DLMO Estimate vs. Sampling Interval",
    x = "Sampling Interval (minutes)",
    y = "Change in DLMO (hours, relative to original)",
    color = "Profile"
  ) +
  theme_minimal()+
  theme(legend.position = "none")


all_results %>%
  filter(resample_interval == 30) %>%
  summarise(
    n = n(),
    all_zero = all(dlmo_difference == 0),
    min_val = min(dlmo_difference),
    max_val = max(dlmo_difference)
  )

library(ggridges)

ggplot(all_results, aes(x = dlmo_difference, y = factor(resample_interval))) +
  geom_density_ridges(
    stat = "binline",
    binwidth = 0.05,
    scale = 1.2,
    draw_baseline = TRUE,
    fill = "skyblue",
    color = "black"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  labs(
    title = "DLMO Change Distributions by Sampling Interval (Histogram Ridges)",
    x = "Change in DLMO (hours)",
    y = "Sampling Interval (minutes)"
  ) +
  theme_minimal()



library(ggridges)

ggplot(all_results, aes(x = dlmo_difference, y = factor(resample_interval), fill = ..x..)) +
  geom_density_ridges_gradient(scale = 1.2, rel_min_height = 0.01) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  labs(
    title = "DLMO Change Distributions by Sampling Interval",
    x = "Change in DLMO (hours)",
    y = "Sampling Interval (minutes)"
  ) +
  theme_minimal()

# mean and SD line plot
all_results %>%
  group_by(resample_interval) %>%
  summarise(
    mean_diff = mean(dlmo_difference),
    sd_diff = sd(dlmo_difference)
  ) %>%
  ggplot(aes(x = resample_interval, y = mean_diff)) +
  geom_line() +
  geom_ribbon(aes(ymin = mean_diff - sd_diff, ymax = mean_diff + sd_diff), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Mean Change in DLMO vs. Sampling Interval",
       x = "Sampling Interval (minutes)", y = "Change in DLMO (hours)")


# prettier
library(ggplot2)
library(dplyr)

# Summarize if not done yet
summary_df <- all_results %>%
  group_by(resample_interval) %>%
  summarise(
    mean_diff = mean(dlmo_difference),
    sd_diff = sd(dlmo_difference),
    n = n(),
    se = sd_diff / sqrt(n),
    .groups = "drop"
  )

# Pretty plot
ggplot(summary_df, aes(x = resample_interval, y = mean_diff)) +
  geom_line(color = "#1f77b4", size = 1) +
  geom_point(size = 2, color = "#1f77b4") +
  geom_ribbon(aes(ymin = mean_diff - se, ymax = mean_diff + se),
              alpha = 0.2, fill = "#1f77b4") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_x_continuous(breaks = summary_df$resample_interval) +
  labs(
    title = "Mean Change in DLMO vs. Sampling Interval",
    x = "Sampling Interval (minutes)",
    y = "Change in DLMO (hours relative to 30-min)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# both SD and SE
library(dplyr)

summary_df <- all_results %>%
  group_by(resample_interval) %>%
  summarise(
    mean_diff = mean(dlmo_difference),
    sd_diff = sd(dlmo_difference),
    n = n(),
    se_diff = sd_diff / sqrt(n),
    .groups = "drop"
  )
library(ggplot2)

ggplot() +
  # SD ribbon (wide, light blue)
  geom_ribbon(data = summary_df,
              aes(x = resample_interval,
                  ymin = mean_diff - sd_diff,
                  ymax = mean_diff + sd_diff),
              fill = "#aec7e8", alpha = 0.3) +

  # SE ribbon (narrow, darker blue)
  geom_ribbon(data = summary_df,
              aes(x = resample_interval,
                  ymin = mean_diff - se_diff,
                  ymax = mean_diff + se_diff),
              fill = "#1f77b4", alpha = 0.3) +

  # Mean line
  geom_line(data = summary_df,
            aes(x = resample_interval, y = mean_diff),
            color = "#1f77b4", size = 1) +

  # Mean points
  geom_point(data = summary_df,
             aes(x = resample_interval, y = mean_diff),
             color = "#1f77b4", size = 2) +

  # Individual data points
  geom_jitter(data = all_results,
              aes(x = resample_interval, y = dlmo_difference),
              width = 1, height = 0, alpha = 0.2, size = 1, color = "black") +

  # Zero reference line
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +

  # Axis & labels
  scale_x_continuous(breaks = summary_df$resample_interval) +
  labs(
    title = "Change in DLMO vs. Sampling Interval",
    subtitle = "Blume et al. results ; Original sampling interval: 30-min",
    x = "Sampling Interval (minutes)",
    y = "Change in DLMO (hours relative to original)"
  ) +

  # Legend-style annotations
  # Horizontal legend annotation at top-right
  annotate("rect", xmin = 65, xmax = 70, ymin = 0.4, ymax = 0.43, fill = "#1f77b4", alpha = 0.3) +  # SE
  annotate("text", x = 71, y = 0.415, label = "Mean ± SE", hjust = 0, size = 4, color = "#1f77b4") +

  annotate("rect", xmin = 65, xmax = 70, ymin = 0.36, ymax = 0.39, fill = "#aec7e8", alpha = 0.3) +  # SD
  annotate("text", x = 71, y = 0.375, label = "± SD", hjust = 0, size = 4, color = "#1f77b4") +


  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

library(dplyr)

ribbons_df <- summary_df %>%
  select(resample_interval, mean_diff, sd_diff, se_diff) %>%
  tidyr::pivot_longer(cols = c(sd_diff, se_diff),
                      names_to = "type", values_to = "spread") %>%
  mutate(
    ymin = mean_diff - spread,
    ymax = mean_diff + spread,
    type = recode(type, sd_diff = "± SD", se_diff = "Mean ± SE")
  )

ggplot() +
  # Ribbons with legend
  geom_ribbon(data = ribbons_df,
              aes(x = resample_interval, ymin = ymin, ymax = ymax, fill = type),
              alpha = 0.3) +

  # Mean line and points
  geom_line(data = summary_df,
            aes(x = resample_interval, y = mean_diff),
            color = "#1f77b4", size = 1) +
  geom_point(data = summary_df,
             aes(x = resample_interval, y = mean_diff),
             color = "#1f77b4", size = 2) +

  # Raw data
  geom_jitter(data = all_results,
              aes(x = resample_interval, y = dlmo_difference),
              width = 1, height = 0, alpha = 0.2, size = 1, color = "black") +

  # Reference line
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +

  scale_x_continuous(breaks = summary_df$resample_interval) +
  scale_fill_manual(values = c("± SD" = "#aec7e8", "Mean ± SE" = "#1f77b4")) +

  labs(
    title = "Change in DLMO vs. Sampling Interval",
    subtitle = "Blume et al. dataset, original sampling interval = 30 min",
    x = "Sampling Interval (minutes)",
    y = "Change in DLMO (decimal hours relative to original)",
    fill = "Shading"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )+
  theme(
    legend.position = c(0.98, 1),       # top right in normalized coordinates
    legend.justification = c(1, 1),        # anchor top-right corner of legend box
    legend.direction = "horizontal",       # horizontal layout
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.background = element_rect(fill = "white", color = NA)
  )




