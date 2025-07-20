# DLMO Monte Carlo Noise Simulation Pipeline
# -------------------------------------------
# Author: [Your Name]
# Date: 2025-07-16
# Purpose: Simulate the impact of time sampling jitter and melatonin measurement noise
#          on Dim Light Melatonin Onset (DLMO) estimates across multiple profiles.

# --- Libraries ----------------------------------------------------------
library(readr)      # fast CSV reading
library(dplyr)      # data manipulation
library(purrr)      # functional mapping
library(fs)         # filesystem operations
library(tibble)     # tidy tibbles
library(dlmoR)      # DLMO estimation functions
library(lubridate)  # date-time handling
library(glue)       # string interpolation
library(future)     # parallel plan
library(furrr)      # parallel map with futures

# --- Utility Functions -------------------------------------------------

# Convert POSIXct times to decimal hours since first observation
posixct_to_decimal <- function(posixct_vec, ref = NULL) {
  if (is.null(ref)) ref <- posixct_vec
  as.numeric(difftime(posixct_vec, min(ref), units = "hours"))
}

# Convert decimal hours back to POSIXct, using the earliest ref_time as origin
decimal_to_posixct <- function(decimal_vec, ref_time) {
  origin_time <- min(ref_time)
  origin_time + dhours(decimal_vec)
}

# Ensure a minimum spacing (in hours) between consecutive time points
enforce_min_time_gap <- function(time_vec, min_gap = 1 / 60) {
  for (i in seq(2, length(time_vec))) {
    if ((time_vec[i] - time_vec[i - 1]) < min_gap) {
      time_vec[i] <- time_vec[i - 1] + min_gap
    }
  }
  time_vec
}

# --- Noise and DLMO Calculation Wrappers --------------------------------

# Safely calculate DLMO without noise; wraps calculate_dlmo() to catch errors
safe_calc_dlmo <- safely(calculate_dlmo)

# Add time jitter and melatonin noise, then estimate DLMO
# Returns a list: $dlmo (the result object) and $perturbed_profile (tibble)
safe_add_noise_and_calc <- function(profile_df, time_sd, mel_sd) {
  safely(function() {
    # 1. Convert to decimal hours
    dh_time <- posixct_to_decimal(profile_df$datetime)

    # 2. Create monotonic spline interpolator of original data
    interp_fun <- splinefun(x = dh_time, y = profile_df$melatonin, method = "monoH.FC")

    # 3. Perturb timepoints with Gaussian noise
    time_pert <- dh_time + rnorm(length(dh_time), 0, time_sd)
    time_sorted <- enforce_min_time_gap(sort(time_pert), min_gap = 1 / 60)

    # 4. Interpolate melatonin at jittered times
    mel_interp <- interp_fun(time_sorted)

    # 5. Add Gaussian noise to melatonin values
    mel_pert <- mel_interp + rnorm(length(mel_interp), 0, mel_sd)

    # 6. Assemble perturbed profile
    perturbed_df <- tibble(
      datetime = decimal_to_posixct(time_sorted, profile_df$datetime),
      melatonin = mel_pert
    )

    # 7. Estimate DLMO on perturbed profile
    list(
      dlmo = calculate_dlmo(perturbed_df),
      perturbed_profile = perturbed_df
    )
  })()
}

# --- Per-Profile Processing ---------------------------------------------

# Reads a CSV, runs baseline + k noisy repeats, and returns a list of results
process_single_profile <- function(profile_path, k, time_sd, mel_sd, output_dir) {
  profile_name <- tools::file_path_sans_ext(basename(profile_path))
  profile_df <- read_csv(profile_path, show_col_types = FALSE)
  message(glue("📄 {profile_name} — baseline + {k} noisy repetitions"))

  result_list <- list()
  flat_results <- list()

  # 1. Baseline (no noise)
  baseline_try <- safe_calc_dlmo(profile_df)
  if (!is.null(baseline_try$result)) {
    res <- baseline_try$result
    row <- tibble(
      profile = profile_name,
      type = "baseline",
      rep = NA_integer_,
      dlmo_decimal = res$ip$inflection_point_fine$x,
      dlmo_result = list(res),
      perturbed_profile = list(NULL),
      time_sd = time_sd,
      mel_sd = mel_sd,
      error_message = NA_character_
    )
  } else {
    row <- tibble(
      profile = profile_name,
      type = "baseline",
      rep = NA_integer_,
      dlmo_decimal = NA_real_,
      dlmo_result = list(NULL),
      perturbed_profile = list(NULL),
      time_sd = time_sd,
      mel_sd = mel_sd,
      error_message = baseline_try$error$message
    )
  }
  result_list[[1]] <- list(type = "baseline", result = row)
  flat_results[[1]] <- row

  # 2. Noisy repetitions
  for (i in seq_len(k)) {
    noise_try <- safe_add_noise_and_calc(profile_df, time_sd, mel_sd)
    if (!is.null(noise_try$result$dlmo)) {
      res <- noise_try$result$dlmo
      perturbed_df <- noise_try$result$perturbed_profile
      row <- tibble(
        profile = profile_name,
        type = "noisy",
        rep = i,
        dlmo_decimal = res$ip$inflection_point_fine$x,
        dlmo_result = list(res),
        perturbed_profile = list(perturbed_df),
        time_sd = time_sd,
        mel_sd = mel_sd,
        error_message = NA_character_
      )
    } else {
      row <- tibble(
        profile = profile_name,
        type = "noisy",
        rep = i,
        dlmo_decimal = NA_real_,
        dlmo_result = list(NULL),
        perturbed_profile = list(NULL),
        time_sd = time_sd,
        mel_sd = mel_sd,
        error_message = noise_try$error$message
      )
    }
    result_list[[length(result_list) + 1]] <- list(type = "noisy", result = row)
    flat_results[[length(flat_results) + 1]] <- row
  }

  # Save per-profile RDS list
  saveRDS(result_list,
          file = file.path(output_dir, paste0(profile_name, "_dlmo_repeats.rds")))

  flat_results
}

# --- Master Pipeline ----------------------------------------------------

# Runs all CSV profiles in a directory, with specified noise levels and repeats
process_all_csv_profiles <- function(
    dir_path = "inst/extdata",
    time_sd = 0.005,
    mel_sd = 0.3,
    k = 25,
    output_dir = "noisy_dlmo_results",
    composite_rds = "all_dlmo_results.rds",
    summary_csv = "all_dlmo_summary.csv"
) {
  dir_create(output_dir)
  profile_files <- list.files(dir_path, pattern = "\\.csv$", full.names = TRUE)
  profile_files <- profile_files[1:2]  # only use first two profiles

  plan(multisession, workers = parallel::detectCores() - 1)

  all_results_flat <- future_map(
    profile_files,
    ~process_single_profile(.x, k, time_sd, mel_sd, output_dir)
  ) %>%
    flatten_df()

  # Save composite results
  saveRDS(all_results_flat, composite_rds)
  write_csv(select(all_results_flat, -dlmo_result, -perturbed_profile), summary_csv)

  message("✅ All results saved:")
  message("  └─ Per-profile RDS → ", output_dir)
  message("  └─ Composite RDS    → ", composite_rds)
  message("  └─ Summary CSV      → ", summary_csv)
}

# --- Simulation Loop over Time Jitter Levels ----------------------------

time_sd_minutes <- c(0, 1, 3, 5, 10)  # jitter levels to test (in minutes)
mel_sd <- 0.3                         # constant melatonin noise (ng/mL)
k <- 2                              # number of noisy repetitions

dir_path <- "inst/extdata"

for (sd_min in time_sd_minutes) {
  time_sd <- sd_min / 60  # convert minutes to hours

  suffix <- glue("{sd_min}min")
  output_dir <- glue("noisy_dlmo_results_{suffix}")
  composite_rds <- glue("all_dlmo_results_{suffix}.rds")
  summary_csv <- glue("all_dlmo_summary_{suffix}.csv")

  message(glue("\n🔁 Running DLMO simulation with time_sd = {sd_min} minutes"))

  process_all_csv_profiles(
    dir_path = dir_path,
    time_sd = time_sd,
    mel_sd = mel_sd,
    k = k,
    output_dir = output_dir,
    composite_rds = composite_rds,
    summary_csv = summary_csv
  )
}

# End of script

# DLMO Noise Impact Plotting Script
# -----------------------------------
# Loads summary CSVs and visualizes the effect of time jitter on DLMO estimation

library(glue)

# ---- Load and Combine Summaries ----------------------------------------

# Specify jitter levels tested (in minutes)
time_sd_minutes <- c(0, 1, 3, 5, 10)
summary_files <- glue("all_dlmo_summary_{time_sd_minutes}min.csv")

# Combine all summaries into one dataframe
all_summaries <- map2_dfr(
  summary_files,
  time_sd_minutes,
  ~ read_csv(.x, show_col_types = FALSE) %>% mutate(time_sd_minutes = .y)
)

# ---- Extract Baseline DLMO for Each Profile ----------------------------

baseline_df <- all_summaries %>%
  filter(type == "baseline") %>%
  select(profile, baseline_dlmo = dlmo_decimal)

# Merge baseline back in for reference
summaries_with_baseline <- all_summaries %>%
  left_join(baseline_df, by = "profile")

# ---- 1. Boxplot of DLMO Estimates by Jitter Level ----------------------

ggplot(
  summaries_with_baseline %>% filter(type == "noisy"),
  aes(x = factor(time_sd_minutes), y = dlmo_decimal)
) +
  geom_boxplot(fill = "skyblue", alpha = 0.6) +
  geom_point(
    data = summaries_with_baseline %>% filter(type == "baseline"),
    aes(x = factor(time_sd_minutes), y = baseline_dlmo),
    color = "red", shape = 18, size = 3
  ) +
  labs(
    title = "Impact of Time Jitter on DLMO Estimates",
    x = "Time jitter (minutes)",
    y = "Estimated DLMO (decimal hours)"
  ) +
  theme_minimal()

# ---- 2. Line Plot of DLMO SD vs. Jitter Level --------------------------

sd_plot_data <- summaries_with_baseline %>%
  filter(type == "noisy", !is.na(dlmo_decimal)) %>%
  group_by(profile, time_sd_minutes) %>%
  summarise(dlmo_sd = sd(dlmo_decimal), .groups = "drop")

ggplot(sd_plot_data, aes(x = time_sd_minutes, y = dlmo_sd, color = profile)) +
  geom_line() +
  geom_point() +
  labs(
    title = "DLMO Estimate Variability vs. Time Jitter",
    x = "Time jitter (minutes)",
    y = "SD of DLMO estimate"
  ) +
  theme_minimal()

# ---- 3. Failure Rate Plot ----------------------------------------------

failure_data <- summaries_with_baseline %>%
  filter(type == "noisy") %>%
  group_by(profile, time_sd_minutes) %>%
  summarise(failure_rate = mean(is.na(dlmo_decimal)), .groups = "drop")

ggplot(failure_data, aes(x = time_sd_minutes, y = failure_rate, color = profile)) +
  geom_line() +
  geom_point() +
  labs(
    title = "DLMO Estimation Failure Rate vs. Time Jitter",
    x = "Time jitter (minutes)",
    y = "Failure rate"
  ) +
  theme_minimal()

# ---- 4. Density Plot of DLMO Distributions -----------------------------

ggplot(
  summaries_with_baseline %>% filter(type == "noisy", !is.na(dlmo_decimal)),
  aes(x = dlmo_decimal)
) +
  geom_density(fill = "lightgray") +
  facet_wrap(~ time_sd_minutes, scales = "free_y") +
  labs(
    title = "Distribution of DLMO Estimates by Jitter Level",
    x = "Estimated DLMO (decimal hours)",
    y = "Density"
  ) +
  theme_minimal()

ggplot(
  summaries_with_baseline %>% filter(type == "noisy", !is.na(dlmo_decimal)),
  aes(x = factor(time_sd_minutes), y = dlmo_decimal)
) +
  geom_violin(fill = "lightblue", color = "black", alpha = 0.6, trim = FALSE) +
  geom_point(
    data = summaries_with_baseline %>% filter(type == "baseline"),
    aes(x = factor(time_sd_minutes), y = baseline_dlmo),
    color = "red", shape = 18, size = 3
  ) +
  labs(
    title = "Distribution of DLMO Estimates by Time Jitter Level",
    x = "Time jitter (minutes)",
    y = "Estimated DLMO (decimal hours)"
  ) +
  theme_minimal()

library(tidyverse)
library(glue)

# --- Load and prepare data (assumes already loaded into summaries_with_baseline) ---

# Filter noisy data only (jitter > 0)
noisy_data <- summaries_with_baseline %>%
  filter(type == "noisy", !is.na(dlmo_decimal))

# Baseline data (jitter = 0)
baseline_points <- summaries_with_baseline %>%
  filter(type == "baseline") %>%
  mutate(time_sd_minutes = 0)

# Combine for unified plotting
plot_data <- bind_rows(noisy_data, baseline_points)

# --- Violin + points plot ---
ggplot(plot_data, aes(x = factor(time_sd_minutes), y = dlmo_decimal)) +
  geom_violin(
    data = noisy_data,
    aes(fill = profile),
    color = "black",
    alpha = 0.3,
    width = 1
  ) +
  geom_jitter(
    data = plot_data,
    aes(color = profile),
    width = 0.15,
    size = 2,
    height = 0
  ) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Impact of Time Jitter on DLMO Estimates",
    x = "Time jitter (minutes)",
    y = "Estimated DLMO (decimal hours)"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

ggplot(
  plot_data %>% filter(type == "noisy", !is.na(dlmo_decimal)),
  aes(x = dlmo_decimal, fill = profile)
) +
  geom_histogram(bins = 50, color = "white", alpha = 0.6) +
  facet_wrap(~ time_sd_minutes, scales = "free_y") +
  labs(
    title = "Distribution of DLMO Estimates by Time Jitter Level",
    x = "Estimated DLMO (decimal hours)",
    y = "Count"
  ) +
  theme_minimal()


# Optional: Save any of the plots using ggsave()
# ggsave("dlmo_jitter_boxplot.png", width = 8, height = 5)

