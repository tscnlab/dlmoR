library(dplyr)
library(dlmoR)
library(purrr)

# --- Determine profiles to process ------------------------------------------
profile_folder <- system.file("extdata/", package = "dlmoR")
csv_files      <- list.files(profile_folder, pattern = "\\.csv$", full.names = TRUE)

# Directory for saving results
out_dir <- "noisy_dlmo_results"

# Identify already processed profiles
result_files  <- list.files(out_dir, pattern = "\\.rds$", full.names = TRUE)
processed_ids  <- tools::file_path_sans_ext(basename(result_files))

# Filter to unprocessed CSVs, naming by profile ID
profiles <- set_names(csv_files,
                      tools::file_path_sans_ext(basename(csv_files)))
profiles <- profiles[!names(profiles) %in% processed_ids]

# --- Noise parameters --------------------------------------------------------
mel_sd_default <- 0.3    # melatonin noise SD
n_rep_default  <- 20     # default repetitions

# --- Helper to enforce minimum time gap --------------------------------------
enforce_min_gap <- function(times_dec, min_gap = 1/60) {
  for (i in seq(2, length(times_dec))) {
    if ((times_dec[i] - times_dec[i - 1]) < min_gap) {
      times_dec[i] <- times_dec[i - 1] + min_gap
    }
  }
  times_dec
}

# --- Function: process profiles ----------------------------------------------
# Loops through each CSV, computes clean + noisy DLMO, and saves results.
process_profiles <- function(profiles,
                             out_dir = "dlmo_results",
                             n_rep   = n_rep_default,
                             mel_sd  = mel_sd_default) {
  dir.create(out_dir, showWarnings = FALSE)

  for (path in profiles) {
    prof_id <- tools::file_path_sans_ext(basename(path))
    message("Processing profile: ", prof_id)

    # Read raw data
    df <- read.csv(path) %>%
      mutate(
        datetime  = as.POSIXct(datetime, tz = "UTC"),
        melatonin = melatonin
      )

    # Compute clean DLMO
    clean_res  <- calculate_dlmo(df)
    clean_time <- clean_res$ip$inflection_point_fine$x

    # Prepare for noise
    dec_t      <- posixct_to_decimal(df$datetime, df$datetime)
    interp_fun <- splinefun(x = dec_t, y = df$melatonin, method = "monoH.FC")
    min_int    <- min(diff(sort(dec_t)))
    time_sd_main <- min_int / 2
    extra_sds    <- c(time_sd_main/2, time_sd_main*2)

    # Helper to run one noise replicate
    run_one <- function(time_sd, mel_sd) {
      t_pert   <- enforce_min_gap(sort(dec_t + rnorm(length(dec_t), 0, time_sd)))
      mel_int  <- interp_fun(t_pert)
      mel_pert <- mel_int + rnorm(length(mel_int), 0, mel_sd)

      sim_df <- tibble(
        datetime  = decimal_to_posixct(t_pert, df$datetime),
        melatonin = mel_pert
      )

      result <- tryCatch(calculate_dlmo(sim_df), error = function(e) NULL)
      xval   <- if (!is.null(result)) result$ip$inflection_point_fine$x else NA_real_
      tibble(
        profile   = prof_id,
        time_sd   = time_sd,
        mel_sd    = mel_sd,
        dlmo_est  = xval
      )
    }

    # Define noise conditions
    conds <- list(
      clean     = list(time_sd = 0,            mel_sd = 0),
      mel_only  = list(time_sd = 0,            mel_sd = mel_sd),
      time_main = list(time_sd = time_sd_main, mel_sd = 0),
      time_ex1  = list(time_sd = extra_sds[1], mel_sd = 0),
      time_ex2  = list(time_sd = extra_sds[2], mel_sd = 0),
      both_axes = list(time_sd = time_sd_main, mel_sd = mel_sd)
    )

    # Run repetitions and collect
    results <- map_df(names(conds), function(cond) {
      p    <- conds[[cond]]
      reps <- if (cond == "clean") 1 else n_rep
      map_dfr(seq_len(reps), ~ run_one(p$time_sd, p$mel_sd) %>% mutate(condition = cond, replicate = .x))
    })

    # Save per-profile
    saveRDS(results, file.path(out_dir, paste0(prof_id, "_dlmo_simple.rds")))
  }
}

# --- Example calls ---------------------------------------------------------
# Run on all unprocessed profiles with 20 reps and default melatonin SD:
# process_profiles(profiles, out_dir = out_dir)
#
# Run on only the first two profiles with 2 replicates each:
 process_profiles(profiles[1:4], out_dir = out_dir, n_rep = 10)
#
# Run on a named subset, e.g. only profiles "subjA" and "subjB":
# selected <- profiles[c("subjA", "subjB")]
# process_profiles(selected, out_dir = out_dir, n_rep = 10, mel_sd = 0.5)

# --- Plotting Results ------------------------------------------------------- -------------------------------------------------------
# library(ggplot2)
# library(dplyr)
#
# all_rds <- list.files(out_dir, pattern = "_dlmo_simple\\.rds$", full.names = TRUE)
# all_res <- purrr::map_df(all_rds, readRDS)
#
# baseline_df <- all_res %>%
#   filter(condition == 'clean') %>%
#   select(profile, baseline = dlmo_est)
#
# plot_df <- all_res %>%
#   filter(condition != 'clean') %>%
#   left_join(baseline_df, by = 'profile') %>%
#   mutate(
#     delta = dlmo_est - baseline,
#     noise_group = case_when(
#       condition == 'mel_only'  ~ 'Melatonin Only',
#       condition == 'both_axes' ~ 'Melatonin + Time',
#       grepl('time', condition) ~ 'Time Only',
#       TRUE                     ~ 'Other'
#     )
#   )
#
# ggplot(plot_df, aes(x = condition, y = delta)) +
#   geom_boxplot() +
#   facet_wrap(~ noise_group, scales = 'free_x') +
#   labs(title = 'Delta DLMO by Noise Condition', x = 'Condition', y = 'Delta (hours)') +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
#

# --- Raincloud Plot (using gghalves) --------------------------------------
# Requires gghalves + ggplot2
# install.packages('gghalves') if needed
# library(gghalves)
#
# # Create raincloud plot: half violin + jittered points + boxplot
# ggplot(plot_df, aes(x = condition, y = delta, fill = noise_group)) +
#   gghalves::geom_half_violin(
#     side = "l",     # left side violin
#     alpha = 0.6,
#     width = 0.6,
#     trim = FALSE
#   ) +
#   geom_boxplot(
#     width = 0.12,
#     outlier.shape = NA,
#     position = position_nudge(x = 0.15)
#   ) +
#   geom_half_point(
#     side = "r",    # points on right side
#     range_scale = 0.4,
#     position = position_nudge(x = 0.15),
#     alpha = 0.4,
#     size = 1
#   ) +
#   coord_flip() +
#   facet_wrap(~ noise_group, scales = 'free_y') +
#   labs(
#     title = 'Raincloud Plot of Delta DLMO by Noise Condition',
#     x = 'Condition',
#     y = 'Delta DLMO (hours)'
#   ) +
#   theme_minimal()

# --- Raincloud Plot (using gghalves) --------------------------------------
# Requires gghalves + ggplot2
# install.packages('gghalves') if needed
# library(gghalves)
#
# plot_df <- all_res %>%
#   filter(condition != "clean") %>%
#   left_join(baseline_df, by = "profile") %>%
#   mutate(
#     delta       = dlmo_est - baseline,
#     noise_group = case_when(
#       condition == "mel_only"  ~ "Melatonin Only",
#       condition == "both_axes" ~ "Melatonin + Time",
#       grepl("time", condition) ~ "Time Only",
#       TRUE                     ~ "Other"
#     ),
#     noise_label = paste0(
#       noise_group,
#       "\n(time SD=", round(time_sd, 2),
#       ", mel SD=",  round(mel_sd, 2), ")"
#     )
#   )
#
# # Vertical raincloud faceted by noise level
# # We facet by noise_label (which includes the SDs)
# ggplot(plot_df, aes(x = condition, y = delta, fill = noise_group)) +
#   gghalves::geom_half_violin(
#     side = 'l',     # left-side violin
#     alpha = 0.6,
#     width = 0.6,
#     trim = FALSE
#   ) +
#   geom_boxplot(
#     width = 0.12,
#     outlier.shape = NA,
#     position = position_nudge(x = 0.15)
#   ) +
#   geom_half_point(
#     side = 'r',    # jittered points on right side
#     range_scale = 0.4,
#     position = position_nudge(x = 0.15),
#     alpha = 0.4,
#     size = 1
#   ) +
#   facet_wrap(~ noise_label, scales = 'free_x') +
#   labs(
#     title = 'Raincloud Plot of Delta DLMO by Noise Condition',
#     x     = 'Condition',
#     y     = 'Delta DLMO (hours)'
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     legend.position = 'none'
#   )

##
# --- Plotting Results: Raincloud Only ---------------------------------------

# --- Plotting Results: Raincloud Only ---------------------------------------

library(ggplot2)
library(dplyr)
library(gghalves)

# Read all saved simple results
all_rds <- list.files(out_dir, pattern = "_dlmo_simple\\.rds$", full.names = TRUE)
all_res <- purrr::map_df(all_rds, readRDS)

# Compute clean baseline per profile
gbaseline <- all_res %>%
  filter(condition == "clean") %>%
  select(profile, baseline = dlmo_est)

# Build plot_df with labels including SD values
plot_df <- all_res %>%
  filter(condition != "clean") %>%
  left_join(gbaseline, by = "profile") %>%
  mutate(
    delta       = dlmo_est - baseline,
    noise_group = case_when(
      condition == "mel_only"  ~ "Melatonin Only",
      condition == "both_axes" ~ "Melatonin + Time",
      grepl("time", condition) ~ "Time Only",
      TRUE                      ~ "Other"
    ),
    noise_label = paste0(
      noise_group,
      "
(time SD = ", round(time_sd, 2),
      ", mel SD = ",  round(mel_sd, 2),
      ")"
    )
  )

# Create a vertical raincloud plot faceted by noise_label
raincloud_plot <- ggplot(plot_df, aes(x = condition, y = delta, fill = noise_group)) +
  gghalves::geom_half_violin(
    side = "l",
    trim = FALSE,
    alpha = 0.6,
    width = 0.6
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    position = position_nudge(x = 0.15)
  ) +
  gghalves::geom_half_point(
    side = "r",
    range_scale = 0.4,
    position = position_nudge(x = 0.15),
    alpha = 0.4,
    size = 1
  ) +
  facet_wrap(~ noise_label, scales = "free_x") +
  labs(
    title = "Raincloud Plot of ΔDLMO by Noise Condition",
    x = "Noise Condition",
    y = "Δ DLMO (hours)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

print(raincloud_plot)

######
# summaries
######

library(dplyr)
library(purrr)

# 1. Read all of your per‑profile result files
all_res <- list.files(out_dir,
                      pattern = "_dlmo_simple\\.rds$",
                      full.names = TRUE) %>%
  map_df(readRDS)

# 2. Peek at the raw tibble
head(all_res)
# or in RStudio use
# View(all_res)

# all_res has columns:
#   profile   – profile ID
#   condition – noise condition (clean, mel_only, time_main, time_ex1, etc.)
#   replicate – replicate number (1 for clean; 1..n_rep for noisy)
#   time_sd, mel_sd, dlmo_est

# 3. Compute “delta” relative to the clean DLMO
baseline <- all_res %>%
  filter(condition == "clean") %>%
  select(profile, baseline = dlmo_est)

full_df <- all_res %>%
  filter(condition != "clean") %>%
  left_join(baseline, by = "profile") %>%
  mutate(delta = dlmo_est - baseline)

# 4. Summarize by profile & condition
summary_df <- full_df %>%
  group_by(profile, condition) %>%
  summarise(
    n            = n(),
    mean_dlmo    = mean(dlmo_est, na.rm = TRUE),
    sd_dlmo      = sd(dlmo_est,   na.rm = TRUE),
    mean_delta   = mean(delta,    na.rm = TRUE),
    sd_delta     = sd(delta,      na.rm = TRUE),
    .groups = "drop"
  )

print(summary_df)
# or View(summary_df)

#######
# deeper summary
#######
library(dplyr)
library(purrr)

# 1. Read all results
all_res <- list.files(out_dir, pattern = "_dlmo_simple\\.rds$", full.names = TRUE) %>%
  map_df(readRDS)

# 2. Extract and join the clean baseline back into every row
baseline <- all_res %>%
  filter(condition == "clean") %>%
  select(profile, baseline = dlmo_est)

full_all <- all_res %>%
  left_join(baseline, by = "profile") %>%
  # delta will be zero for the clean run
  mutate(delta = dlmo_est - baseline)

# Now `full_all` has one row per profile × condition × replicate, with columns:
#   profile, condition, replicate, time_sd, mel_sd, dlmo_est, baseline, delta

# Inspect the individual-level table
print(full_all)
# Or in RStudio: View(full_all)

# 3. Summary by profile & condition
summary_df <- full_all %>%
  group_by(profile, condition) %>%
  summarise(
    n_reps      = n(),
    mean_dlmo   = mean(dlmo_est, na.rm = TRUE),
    sd_dlmo     = sd(dlmo_est,   na.rm = TRUE),
    mean_delta  = mean(delta,    na.rm = TRUE),
    sd_delta    = sd(delta,      na.rm = TRUE),
    .groups     = "drop"
  )

print(summary_df)
# Or: View(summary_df)


