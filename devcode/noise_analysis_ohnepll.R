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
processed_ids <- tools::file_path_sans_ext(basename(result_files))

# Filter to unprocessed CSVs, naming by profile ID
profiles <- set_names(csv_files,
                      tools::file_path_sans_ext(basename(csv_files)))
profiles <- profiles[!names(profiles) %in% processed_ids]

# --- Noise parameters --------------------------------------------------------
cv_intra      <- 7.9      # intra‑assay CV in %
n_rep_default <- 20       # default repetitions

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
                             out_dir = "noisy_dlmo_results",
                             n_rep   = n_rep_default,
                             cv      = cv_intra) {
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

    # Compute clean DLMO for baseline
    clean_res <- calculate_dlmo(df)
    clean_time <- clean_res$ip$inflection_point_fine$x

    # Prepare for interpolation
    dec_t      <- posixct_to_decimal(df$datetime, df$datetime)
    interp_fun <- splinefun(x = dec_t, y = df$melatonin, method = "monoH.FC")

    # Define your fixed time SDs (in hours)
    time_sd_5  <-  5 / 60
    time_sd_10 <- 10 / 60
    time_sd_20 <- 20 / 60

    # Single-run simulator
    run_one <- function(time_sd, add_mel_noise = FALSE) {
      # 1) jitter time
      t_pert   <- enforce_min_gap(sort(dec_t + rnorm(length(dec_t), 0, time_sd)))
      # 2) interpolate the true curve at the jittered times
      mel_int  <- interp_fun(t_pert)
      # 3) optionally add melatonin noise from CV
      # 3) optionally add melatonin noise from CV, and compute SD stats
      if (add_mel_noise) {
        sd_vec       <- (cv / 100) * mel_int
        mel_pert     <- mel_int + rnorm(length(mel_int), 0, sd_vec)
        mean_mel_sd  <- mean(sd_vec)
        min_mel_sd   <- min(sd_vec)
        max_mel_sd   <- max(sd_vec)
      } else {
        mel_pert     <- mel_int
        mean_mel_sd  <- 0
        min_mel_sd   <- 0
        max_mel_sd   <- 0
      }

      sim_df <- tibble(
        datetime  = decimal_to_posixct(t_pert, df$datetime),
        melatonin = mel_pert
      )

      # compute DLMO, catch errors
      result <- tryCatch(calculate_dlmo(sim_df), error = function(e) NULL)
      xval   <- if (!is.null(result)) result$ip$inflection_point_fine$x else NA_real_

      tibble(
        profile     = prof_id,
        time_sd     = time_sd,
        mean_mel_sd = mean_mel_sd,
        min_mel_sd  = min_mel_sd,
        max_mel_sd  = max_mel_sd,
        dlmo_est    = xval
      )
    }

    # Set up all conditions
    conds <- list(
      clean         = list(time_sd = 0,           add_mel = FALSE),
      mel_only      = list(time_sd = 0,           add_mel = TRUE),
      time_5min     = list(time_sd = time_sd_5,   add_mel = FALSE),
      time_10min    = list(time_sd = time_sd_10,  add_mel = FALSE),
      time_20min    = list(time_sd = time_sd_20,  add_mel = FALSE),
      both_axes_10  = list(time_sd = time_sd_10,  add_mel = TRUE)
    )

    # Run sims
    results <- map_df(names(conds), function(cond) {
      p    <- conds[[cond]]
      reps <- if (cond == "clean") 1 else n_rep
      map_dfr(seq_len(reps), function(i) {
        run_one(p$time_sd, p$add_mel) %>%
          mutate(condition = cond,
                 replicate = i)
      })
    })

    # Save per-profile
    saveRDS(results,
            file.path(out_dir, paste0(prof_id, "_dlmo_time_mel_stats.rds")))
  }
}


# --- Example call: run on all unprocessed profiles --------------------------
#process_profiles(profiles, out_dir = out_dir, n_rep = 20)

# Run on only the first two profiles with 2 replicates each:
process_profiles(profiles[1:3], out_dir = out_dir, n_rep = 5)

##############
# prepare summariers
##############


library(dplyr)
library(purrr)
library(ggplot2)
library(gghalves)

# 1. Read in your per‐profile stats files
all_rds <- list.files(out_dir,
                      pattern = "_dlmo_time_mel_stats\\.rds$",
                      full.names = TRUE)
all_res <- map_df(all_rds, readRDS)

# Quick check
print(names(all_res))
# should include: profile, condition, replicate, time_sd,
#                 mean_mel_sd, min_mel_sd, max_mel_sd, dlmo_est

# 2. Extract baseline (clean) per profile
baseline_df <- all_res %>%
  filter(condition == "clean") %>%
  select(profile, baseline = dlmo_est)

# 3. Join into full table and compute delta
full_df <- all_res %>%
  left_join(baseline_df, by = "profile") %>%
  mutate(delta = dlmo_est - baseline)

# 4. Summary table (includes the clean run as one of the reps)
summary_df <- full_df %>%
  group_by(profile, condition) %>%
  summarise(
    n_reps       = n(),
    mean_dlmo    = mean(dlmo_est,   na.rm = TRUE),
    sd_dlmo      = sd(dlmo_est,     na.rm = TRUE),
    mean_delta   = mean(delta,      na.rm = TRUE),
    sd_delta     = sd(delta,        na.rm = TRUE),
    mean_mel_sd  = mean(mean_mel_sd, na.rm = TRUE),
    min_mel_sd   = min(min_mel_sd,  na.rm = TRUE),
    max_mel_sd   = max(max_mel_sd,  na.rm = TRUE),
    mean_time_sd = mean(time_sd,     na.rm = TRUE),
    .groups      = "drop"
  )

print(summary_df)
View(summary_df) #in RStudio for a spreadsheet‐style look

####################
# Raincloud of delta DLMO for each condition
####################
library(dplyr)
library(purrr)
library(ggplot2)
library(gghalves)

# 1) Read in all of your *_dlmo_time_mel_stats.rds files
all_res <- list.files(out_dir,
                      pattern = "_dlmo_time_mel_stats\\.rds$",
                      full.names = TRUE) %>%
  map_df(readRDS)

# 2) Pull out the clean (baseline) DLMO per profile
baseline <- all_res %>%
  filter(condition == "clean") %>%
  select(profile, clean_dlmo = dlmo_est)

# 3) Join & compute delta for *every* row (including clean)
full_df <- all_res %>%
  left_join(baseline, by = "profile") %>%
  mutate(
    delta = dlmo_est - clean_dlmo
  )

# 4) Select exactly the six conditions you care about
wanted <- c("clean","mel_only","time_5min","time_10min","time_20min","both_axes_10")
labels <- c(
  "Clean",
  "Melatonin Only",
  "Time Only\n(5 min)",
  "Time Only\n(10 min)",
  "Time Only\n(20 min)",
  "Mel + Time\n(10 min)"
)

plot_df <- full_df %>%
  filter(condition %in% wanted) %>%
  mutate(
    condition = factor(condition, levels = wanted, labels = labels)
  )

# 5) Split out clean vs noisy
df_clean    <- filter(plot_df, condition == "Clean")
df_no_clean <- filter(plot_df, condition != "Clean")

# 6) Plot—all deltas are pre‑computed, so clean is exactly zero
ggplot() +
  # Raincloud for noisy conditions
  geom_half_violin(
    data = df_no_clean,
    aes(x = condition, y = delta, fill = condition),
    side  = "l", trim = FALSE, alpha = 0.6, width = 0.6
  ) +
  geom_boxplot(
    data = df_no_clean,
    aes(x = condition, y = delta),
    width = 0.10, outlier.shape = NA,
    position = position_nudge(x = 0.1), alpha = 0.8
  ) +
  geom_half_point(
    data        = df_no_clean,
    aes(x = condition, y = delta),
    side        = "r", range_scale = 0.4,
    position    = position_nudge(x = 0.1),
    alpha       = 0.4, size = 1
  ) +
  # Plot the clean runs as points—delta is already zero
  # geom_point(
  #   data = df_clean,
  #   aes(x = condition, y = delta),
  #   color = "black", fill = "black",
  #   shape = 21, size = 3
  # ) +
  labs(
    title = "Raincloud Plot of ΔDLMO by Noise Condition",
    x     = NULL,
    y     = "Δ DLMO (hours)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x       = element_text(angle = 45, hjust = 1),
    legend.position   = "none",
    panel.grid.major.x = element_blank()
  )
#
# #########
# # raincloud of DLMO estimates for each condition
# #########
#
# # library(dplyr)
# # library(purrr)
# # library(ggplot2)
# # library(gghalves)
# #
# # # 1. Read in all of your *_dlmo_time_mel_stats.rds files
# # all_res <- list.files(out_dir,
# #                       pattern = "_dlmo_time_mel_stats\\.rds$",
# #                       full.names = TRUE) %>%
# #   map_df(readRDS)
# #
# # # 2. Compute baseline and delta
# # baseline <- all_res %>%
# #   filter(condition == "clean") %>%
# #   select(profile, baseline = dlmo_est)
# #
# # full_df <- all_res %>%
# #   left_join(baseline, by = "profile") %>%
# #   mutate(delta = dlmo_est - baseline)
# #
# # # 3. Keep only your six conditions and relabel them
# # wanted <- c("clean","mel_only","time_5min","time_10min","time_20min","both_axes_10")
# # labels <- c(
# #   "Clean",
# #   "Melatonin Only",
# #   "Time Only\n(5 min)",
# #   "Time Only\n(10 min)",
# #   "Time Only\n(20 min)",
# #   "Mel + Time\n(10 min)"
# # )
# #
# # plot_df <- full_df %>%
# #   filter(condition %in% wanted) %>%
# #   mutate(
# #     condition = factor(condition, levels = wanted, labels = labels)
# #   )
# #
# # # 4. Draw a pooled raincloud (half‐violin + points + box) per condition
# # ggplot(plot_df, aes(x = condition, y = delta, fill = condition)) +
# #   gghalves::geom_half_violin(
# #     side = "l",
# #     trim = FALSE,
# #     alpha = 0.6,
# #     width = 0.6
# #   ) +
# #   geom_boxplot(
# #     width = 0.1,
# #     outlier.shape = NA,
# #     position = position_nudge(x = 0.1),
# #     alpha = 0.8
# #   ) +
# #   gghalves::geom_half_point(
# #     side = "r",
# #     range_scale = 0.4,
# #     position = position_nudge(x = 0.1),
# #     alpha = 0.4,
# #     size = 1
# #   ) +
# #   labs(
# #     title = "Raincloud Plot of ΔDLMO by Noise Condition",
# #     x     = NULL,
# #     y     = "Δ DLMO (hours)"
# #   ) +
# #   theme_minimal() +
# #   theme(
# #     axis.text.x    = element_text(angle = 45, hjust = 1),
# #     legend.position = "none",
# #     panel.grid.major.x = element_blank()
# #   )
#
#
#
# library(dplyr)
# library(purrr)
# library(ggplot2)
#
# # 1) Read in all per‐profile stats
# all_res <- list.files(out_dir,
#                       pattern = "_dlmo_time_mel_stats\\.rds$",
#                       full.names = TRUE) %>%
#   map_df(readRDS)
#
# # 2) For each profile, compute delta = dlmo_est - clean_dlmo
# delta_df <- all_res %>%
#   group_by(profile) %>%
#   mutate(
#     clean_dlmo = dlmo_est[condition == "clean"],
#     delta      = dlmo_est - clean_dlmo
#   ) %>%
#   ungroup()
#
# # 3) Keep just the six conditions you care about, with nice labels
# wanted <- c("clean","mel_only","time_5min","time_10min","time_20min","both_axes_10")
# labels <- c(
#   "Clean",
#   "Melatonin Only",
#   "Time Only\n(5 min)",
#   "Time Only\n(10 min)",
#   "Time Only\n(20 min)",
#   "Mel + Time\n(10 min)"
# )
#
# plot_df <- delta_df %>%
#   filter(condition %in% wanted) %>%
#   mutate(
#     condition = factor(condition, levels = wanted, labels = labels)
#   )
#
# # 4) Violin‐plot of ΔDLMO by condition (all profiles & reps pooled)
# ggplot(plot_df, aes(x = condition, y = delta, fill = condition)) +
#   geom_violin(trim = FALSE, alpha = 0.6) +
#   geom_boxplot(width = 0.1, outlier.shape = NA, position = position_nudge(x = 0)) +
#   stat_summary(fun = median, geom = "point", size = 2, color = "black") +
#   labs(
#     title = "Violin Plot of ΔDLMO by Noise Condition",
#     x     = NULL,
#     y     = expression(Delta~DLMO~"(hours)")
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.x        = element_text(angle = 45, hjust = 1),
#     legend.position    = "none",
#     panel.grid.major.x = element_blank()
#   )
#
#
#
# # 5. Raincloud plot
# plot_df <- full_df %>%
#   filter(condition != "clean") %>%
#   mutate(
#     noise_group = case_when(
#       condition == "mel_only"     ~ "Melatonin Only",
#       condition == "both_axes_10" ~ "Melatonin + Time (10 min)",
#       grepl("time_5min", condition)  ~ "Time Only (5 min)",
#       grepl("time_10min",condition)  ~ "Time Only (10 min)",
#       grepl("time_20min",condition)  ~ "Time Only (20 min)",
#       TRUE                         ~ condition
#     ),
#     noise_label = paste0(
#       noise_group,
#       "\n(time SD = ", round(time_sd, 2),
#       " h, mel SD = ", round(mean_mel_sd, 2), " pg/mL)"
#     )
#   )
#
# raincloud_plot <- ggplot(plot_df, aes(x = condition, y = delta, fill = noise_group)) +
#   gghalves::geom_half_violin(side = "l", trim = FALSE, alpha = 0.6, width = 0.6) +
#   geom_boxplot(width = 0.12, outlier.shape = NA, position = position_nudge(x = 0.15)) +
#   gghalves::geom_half_point(side = "r", range_scale = 0.4,
#                             position = position_nudge(x = 0.15),
#                             alpha = 0.4, size = 1) +
#   facet_wrap(~ noise_label, scales = "free_x") +
#   labs(
#     title = "Raincloud Plot of ΔDLMO by Noise Condition",
#     x     = "Condition",
#     y     = "Δ DLMO (hours)"
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.x    = element_text(angle = 45, hjust = 1),
#     legend.position = "none"
#   )
#
# print(raincloud_plot)
#
#
# library(dplyr)
# library(purrr)
# library(ggplot2)
#
# # 1) Read in all of your *_dlmo_time_mel_stats.rds files
# all_res <- list.files(out_dir,
#                       pattern = "_dlmo_time_mel_stats\\.rds$",
#                       full.names = TRUE) %>%
#   map_df(readRDS)
#
# # 2) Get the “clean” (no‐noise) baseline per profile
# baseline <- all_res %>%
#   filter(condition == "clean") %>%
#   select(profile, baseline = dlmo_est)
#
# # 3) Join and compute delta for *every* row
# full_df <- all_res %>%
#   left_join(baseline, by = "profile") %>%
#   mutate(delta = dlmo_est - baseline)
#
# # 4) If you want *exactly* one violin per the 6 conditions you defined:
# #    clean, mel_only, time_5min, time_10min, time_20min, both_axes_10
# wanted <- c("clean","mel_only","time_5min","time_10min","time_20min","both_axes_10")
# plot_df <- full_df %>%
#   filter(condition %in% wanted)
#
# # 5) (Optional) give nicer labels and enforce order
# plot_df <- plot_df %>%
#   mutate(
#     condition = factor(condition, levels = wanted,
#                        labels = c(
#                          "Clean",
#                          "Melatonin Only",
#                          "Time Only\n(5 min)",
#                          "Time Only\n(10 min)",
#                          "Time Only\n(20 min)",
#                          "Mel + Time\n(10 min)"
#                        )
#     )
#   )
#
# # 6) Draw your pooled violin plot of ΔDLMO
# ggplot(plot_df, aes(x = condition, y = delta)) +
#   geom_violin(fill = "#A6CEE3", color = NA, alpha = 0.8) +
#   stat_summary(fun = median, geom = "point", size = 2, color = "darkred") +
#   labs(
#     title = "Pooled ΔDLMO by Noise Condition (all profiles)",
#     x     = NULL,
#     y     = "Δ DLMO (hours)"
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     panel.grid.major.x = element_blank()
#   )
# library(dplyr)
# library(purrr)
# library(ggplot2)
# library(gghalves)
#
# # read in your stats
# all_res <- list.files(out_dir,
#                       pattern = "_dlmo_time_mel_stats\\.rds$",
#                       full.names = TRUE) %>%
#   map_df(readRDS)
#
# # keep & label the 6 conditions
# wanted <- c("clean","mel_only","time_5min","time_10min","time_20min","both_axes_10")
# labels <- c("Clean","Melatonin Only","Time Only\n(5 min)","Time Only\n(10 min)",
#             "Time Only\n(20 min)","Mel + Time\n(10 min)")
#
# plot_df <- all_res %>%
#   filter(condition %in% wanted) %>%
#   mutate(condition = factor(condition, levels = wanted, labels = labels))
#
# df_no_clean <- filter(plot_df, condition != "Clean")
# df_clean    <- filter(plot_df, condition == "Clean")
#
# ggplot() +
#   # half‐violins for noisy conditions
#   geom_half_violin(
#     data = df_no_clean,
#     aes(x = condition, y = dlmo_est, fill = condition),
#     side  = "l", trim = FALSE, alpha = 0.6, width = 0.6
#   ) +
#   # boxplots
#   geom_boxplot(
#     data = df_no_clean,
#     aes(x = condition, y = dlmo_est),
#     width = 0.10, outlier.shape = NA,
#     position = position_nudge(x = 0.1), alpha = 0.8
#   ) +
#   # half‐points for noisy: no built‐in jitter, but add x‐only jitter
#   gghalves::geom_half_point(
#     data        = df_no_clean,
#     aes(x = condition, y = dlmo_est),
#     side        = "r",
#     range_scale = 0,  # turn off all internal jitter
#     position    = position_jitter(width = 0.05, height = 0),
#     alpha       = 0.4,
#     size        = 1
#   ) +
#   # points for Clean: same x‐jitter, y exactly dlmo_est
#   geom_point(
#     data     = df_clean,
#     aes(x = condition, y = dlmo_est),
#     position = position_jitter(width = 0.05, height = 0),
#     alpha    = 0.4,
#     size     = 1
#   ) +
#   labs(
#     title = "Raincloud Plot of DLMO Estimates by Noise Condition",
#     x     = NULL,
#     y     = "DLMO Estimate (hours)"
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.x        = element_text(angle = 45, hjust = 1),
#     legend.position    = "none",
#     panel.grid.major.x = element_blank()
#   )
#
#
# #### use this DLMO violins
# library(dplyr)
# library(purrr)
# library(ggplot2)
# library(gghalves)
#
# # 1. Read in your per‐profile stats files
# all_res <- list.files(out_dir,
#                       pattern = "_dlmo_time_mel_stats\\.rds$",
#                       full.names = TRUE) %>%
#   map_df(readRDS)
#
# # 2. Select the six conditions and give them nice labels
# wanted <- c("clean","mel_only","time_5min","time_10min","time_20min","both_axes_10")
# labels <- c(
#   "Clean",
#   "Melatonin Only",
#   "Time Only\n(5 min)",
#   "Time Only\n(10 min)",
#   "Time Only\n(20 min)",
#   "Mel + Time\n(10 min)"
# )
# ggplot(plot_df, aes(x = condition, y = dlmo_est, fill = condition)) +
#   # half‑violins for all conditions
#   gghalves::geom_half_violin(
#     side  = "l",
#     trim  = FALSE,
#     alpha = 0.6,
#     width = 0.6
#   ) +
#   # boxplots
#   geom_boxplot(
#     width         = 0.10,
#     outlier.shape = NA,
#     position      = position_nudge(x = 0.1),
#     alpha         = 0.8
#   ) +
#   # half‑points
#   gghalves::geom_half_point(
#     side        = "r",
#     range_scale = 0,
#     position    = position_nudge(x = 0.1),
#     alpha       = 0.4,
#     size        = 1
#   ) +
#   labs(
#     title = "Raincloud Plot of DLMO Estimates by Noise Condition",
#     x     = NULL,
#     y     = "DLMO Estimate (hours)"
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.x        = element_text(angle = 45, hjust = 1),
#     legend.position    = "none",
#     panel.grid.major.x = element_blank()
#   )
#
#
#
# ## correct delta plot
# library(dplyr)
# library(purrr)
# library(ggplot2)
#
# # 1) Read in all per‐profile stats
# all_res <- list.files(out_dir,
#                       pattern = "_dlmo_time_mel_stats\\.rds$",
#                       full.names = TRUE) %>%
#   map_df(readRDS)
#
# # 2) For each profile, compute delta = dlmo_est - clean_dlmo
# delta_df <- all_res %>%
#   group_by(profile) %>%
#   mutate(
#     clean_dlmo = dlmo_est[condition == "clean"],
#     delta      = dlmo_est - clean_dlmo
#   ) %>%
#   ungroup()
#
# # 3) Keep just the six conditions you care about, with nice labels
# wanted <- c("clean","mel_only","time_5min","time_10min","time_20min","both_axes_10")
# labels <- c(
#   "Clean",
#   "Melatonin Only",
#   "Time Only\n(5 min)",
#   "Time Only\n(10 min)",
#   "Time Only\n(20 min)",
#   "Mel + Time\n(10 min)"
# )
#
# plot_df <- delta_df %>%
#   filter(condition %in% wanted) %>%
#   mutate(
#     condition = factor(condition, levels = wanted, labels = labels)
#   )
#
# # 4) Violin‐plot of ΔDLMO by condition (all profiles & reps pooled)
# ggplot(plot_df, aes(x = condition, y = delta, fill = condition)) +
#   geom_violin(trim = FALSE, alpha = 0.6) +
#   geom_boxplot(width = 0.1, outlier.shape = NA, position = position_nudge(x = 0)) +
#   stat_summary(fun = median, geom = "point", size = 2, color = "black") +
#   labs(
#     title = "Violin Plot of ΔDLMO by Noise Condition",
#     x     = NULL,
#     y     = expression(Delta~DLMO~"(hours)")
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.x        = element_text(angle = 45, hjust = 1),
#     legend.position    = "none",
#     panel.grid.major.x = element_blank()
#   )
#
#
#
# # --- Plotting Results ------------------------------------------------------- -------------------------------------------------------
# # library(ggplot2)
# # library(dplyr)
# #
# # all_rds <- list.files(out_dir, pattern = "_dlmo_simple\\.rds$", full.names = TRUE)
# # all_res <- purrr::map_df(all_rds, readRDS)
# #
# # baseline_df <- all_res %>%
# #   filter(condition == 'clean') %>%
# #   select(profile, baseline = dlmo_est)
# #
# # plot_df <- all_res %>%
# #   filter(condition != 'clean') %>%
# #   left_join(baseline_df, by = 'profile') %>%
# #   mutate(
# #     delta = dlmo_est - baseline,
# #     noise_group = case_when(
# #       condition == 'mel_only'  ~ 'Melatonin Only',
# #       condition == 'both_axes' ~ 'Melatonin + Time',
# #       grepl('time', condition) ~ 'Time Only',
# #       TRUE                     ~ 'Other'
# #     )
# #   )
# #
# # ggplot(plot_df, aes(x = condition, y = delta)) +
# #   geom_boxplot() +
# #   facet_wrap(~ noise_group, scales = 'free_x') +
# #   labs(title = 'Delta DLMO by Noise Condition', x = 'Condition', y = 'Delta (hours)') +
# #   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# #
#
# # --- Raincloud Plot (using gghalves) --------------------------------------
# # Requires gghalves + ggplot2
# # install.packages('gghalves') if needed
# # library(gghalves)
# #
# # # Create raincloud plot: half violin + jittered points + boxplot
# # ggplot(plot_df, aes(x = condition, y = delta, fill = noise_group)) +
# #   gghalves::geom_half_violin(
# #     side = "l",     # left side violin
# #     alpha = 0.6,
# #     width = 0.6,
# #     trim = FALSE
# #   ) +
# #   geom_boxplot(
# #     width = 0.12,
# #     outlier.shape = NA,
# #     position = position_nudge(x = 0.15)
# #   ) +
# #   geom_half_point(
# #     side = "r",    # points on right side
# #     range_scale = 0.4,
# #     position = position_nudge(x = 0.15),
# #     alpha = 0.4,
# #     size = 1
# #   ) +
# #   coord_flip() +
# #   facet_wrap(~ noise_group, scales = 'free_y') +
# #   labs(
# #     title = 'Raincloud Plot of Delta DLMO by Noise Condition',
# #     x = 'Condition',
# #     y = 'Delta DLMO (hours)'
# #   ) +
# #   theme_minimal()
#
# # --- Raincloud Plot (using gghalves) --------------------------------------
# # Requires gghalves + ggplot2
# # install.packages('gghalves') if needed
# # library(gghalves)
# #
# # plot_df <- all_res %>%
# #   filter(condition != "clean") %>%
# #   left_join(baseline_df, by = "profile") %>%
# #   mutate(
# #     delta       = dlmo_est - baseline,
# #     noise_group = case_when(
# #       condition == "mel_only"  ~ "Melatonin Only",
# #       condition == "both_axes" ~ "Melatonin + Time",
# #       grepl("time", condition) ~ "Time Only",
# #       TRUE                     ~ "Other"
# #     ),
# #     noise_label = paste0(
# #       noise_group,
# #       "\n(time SD=", round(time_sd, 2),
# #       ", mel SD=",  round(mel_sd, 2), ")"
# #     )
# #   )
# #
# # # Vertical raincloud faceted by noise level
# # # We facet by noise_label (which includes the SDs)
# # ggplot(plot_df, aes(x = condition, y = delta, fill = noise_group)) +
# #   gghalves::geom_half_violin(
# #     side = 'l',     # left-side violin
# #     alpha = 0.6,
# #     width = 0.6,
# #     trim = FALSE
# #   ) +
# #   geom_boxplot(
# #     width = 0.12,
# #     outlier.shape = NA,
# #     position = position_nudge(x = 0.15)
# #   ) +
# #   geom_half_point(
# #     side = 'r',    # jittered points on right side
# #     range_scale = 0.4,
# #     position = position_nudge(x = 0.15),
# #     alpha = 0.4,
# #     size = 1
# #   ) +
# #   facet_wrap(~ noise_label, scales = 'free_x') +
# #   labs(
# #     title = 'Raincloud Plot of Delta DLMO by Noise Condition',
# #     x     = 'Condition',
# #     y     = 'Delta DLMO (hours)'
# #   ) +
# #   theme_minimal() +
# #   theme(
# #     axis.text.x = element_text(angle = 45, hjust = 1),
# #     legend.position = 'none'
# #   )
#
# ##
# # --- Plotting Results: Raincloud Only ---------------------------------------
#
# # --- Plotting Results: Raincloud Only ---------------------------------------
#
# library(ggplot2)
# library(dplyr)
# library(gghalves)
#
# # Read all saved simple results
# # read in the *time_mel_stats* files, not the *_simple.rds
# all_rds <- list.files(out_dir,
#                       pattern = "_dlmo_time_mel_stats\\.rds$",
#                       full.names = TRUE)
#
# all_res <- purrr::map_df(all_rds, readRDS)
#
# # now inspect what you got:
# names(all_res)
#
#
# # Compute clean baseline per profile
# gbaseline <- all_res %>%
#   filter(condition == "clean") %>%
#   select(profile, baseline = dlmo_est)
#
# # Build plot_df with labels including SD values
# plot_df <- all_res %>%
#   filter(condition != "clean") %>%
#   left_join(gbaseline, by = "profile") %>%
#   mutate(
#     delta       = dlmo_est - baseline,
#     noise_group = case_when(
#       condition == "mel_only"  ~ "Melatonin Only",
#       condition == "both_axes" ~ "Melatonin + Time",
#       grepl("time", condition) ~ "Time Only",
#       TRUE                      ~ "Other"
#     ),
#     noise_label = paste0(
#       noise_group,
#       "
# (time SD = ", round(time_sd, 2),
#       ", mel SD = ",  round(mel_sd, 2),
#       ")"
#     )
#   )
#
# # Create a vertical raincloud plot faceted by noise_label
# raincloud_plot <- ggplot(plot_df, aes(x = condition, y = delta, fill = noise_group)) +
#   gghalves::geom_half_violin(
#     side = "l",
#     trim = FALSE,
#     alpha = 0.6,
#     width = 0.6
#   ) +
#   geom_boxplot(
#     width = 0.12,
#     outlier.shape = NA,
#     position = position_nudge(x = 0.15)
#   ) +
#   gghalves::geom_half_point(
#     side = "r",
#     range_scale = 0.4,
#     position = position_nudge(x = 0.15),
#     alpha = 0.4,
#     size = 1
#   ) +
#   facet_wrap(~ noise_label, scales = "free_x") +
#   labs(
#     title = "Raincloud Plot of ΔDLMO by Noise Condition",
#     x = "Noise Condition",
#     y = "Δ DLMO (hours)"
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     legend.position = "none"
#   )
#
# print(raincloud_plot)
#
# ######
# # summaries
# ######
#
# library(dplyr)
# library(purrr)
#
# # 1. Read all of your per‑profile result files
# all_res <- list.files(out_dir,
#                       pattern = "_dlmo_simple\\.rds$",
#                       full.names = TRUE) %>%
#   map_df(readRDS)
#
# # 2. Peek at the raw tibble
# head(all_res)
# # or in RStudio use
# # View(all_res)
#
# # all_res has columns:
# #   profile   – profile ID
# #   condition – noise condition (clean, mel_only, time_main, time_ex1, etc.)
# #   replicate – replicate number (1 for clean; 1..n_rep for noisy)
# #   time_sd, mel_sd, dlmo_est
#
# # 3. Compute “delta” relative to the clean DLMO
# baseline <- all_res %>%
#   filter(condition == "clean") %>%
#   select(profile, baseline = dlmo_est)
#
# full_df <- all_res %>%
#   filter(condition != "clean") %>%
#   left_join(baseline, by = "profile") %>%
#   mutate(delta = dlmo_est - baseline)
#
# # 4. Summarize by profile & condition
# summary_df <- full_df %>%
#   group_by(profile, condition) %>%
#   summarise(
#     n            = n(),
#     mean_dlmo    = mean(dlmo_est, na.rm = TRUE),
#     sd_dlmo      = sd(dlmo_est,   na.rm = TRUE),
#     mean_delta   = mean(delta,    na.rm = TRUE),
#     sd_delta     = sd(delta,      na.rm = TRUE),
#     .groups = "drop"
#   )
#
# print(summary_df)
# # or View(summary_df)
#
# #######
# # deeper summary
# #######
# library(dplyr)
# library(purrr)
#
# # 1. Read all results
# all_res <- list.files(out_dir, pattern = "_dlmo_simple\\.rds$", full.names = TRUE) %>%
#   map_df(readRDS)
#
# # 2. Extract and join the clean baseline back into every row
# baseline <- all_res %>%
#   filter(condition == "clean") %>%
#   select(profile, baseline = dlmo_est)
#
# full_all <- all_res %>%
#   left_join(baseline, by = "profile") %>%
#   # delta will be zero for the clean run
#   mutate(delta = dlmo_est - baseline)
#
# # Now `full_all` has one row per profile × condition × replicate, with columns:
# #   profile, condition, replicate, time_sd, mel_sd, dlmo_est, baseline, delta
#
# # Inspect the individual-level table
# print(full_all)
# # Or in RStudio: View(full_all)
#
# # 3. Summary by profile & condition
# summary_df <- full_all %>%
#   group_by(profile, condition) %>%
#   summarise(
#     n_reps      = n(),
#     mean_dlmo   = mean(dlmo_est, na.rm = TRUE),
#     sd_dlmo     = sd(dlmo_est,   na.rm = TRUE),
#     mean_delta  = mean(delta,    na.rm = TRUE),
#     sd_delta    = sd(delta,      na.rm = TRUE),
#     .groups     = "drop"
#   )
#
# print(summary_df)
# # Or: View(summary_df)
#
#
