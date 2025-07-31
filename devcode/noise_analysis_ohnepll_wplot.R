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
# all_rds <- list.files(out_dir,
#                       pattern = "_dlmo_time_mel_stats\\.rds$",
#                       full.names = TRUE)
# all_res <- map_df(all_rds, readRDS)
# 1) Read in all files that start with "civibe_melatonin_"
out_dir = "~/Documents/Projects/DLMO/dlmoRpaperresults/Civibe/noise_results_skinny"
out_dir = "~/Documents/Projects/DLMO/dlmoRpaperresults/Blume/noise_results_skinny"

# 1) List all .rds files that start with "civibe_melatonin_"
all_files <- list.files(out_dir,
                        # pattern = "^civibe_melatonin_.*\\.rds$",
                        pattern = "\\.rds$",
                        full.names = TRUE)

all_res <- list.files(out_dir,
                      pattern = "\\.rds",
                      full.names = TRUE) %>%
  map_df(readRDS)

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
# View(full_df)
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

#print(summary_df)
#View(summary_df) #in RStudio for a spreadsheet‐style look

####################
# Raincloud of delta DLMO for each condition
####################
library(dplyr)
library(purrr)
library(ggplot2)
library(gghalves)

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




########
library(dplyr)
library(ggplot2)
library(gghalves)

# Map subgroup to manual numeric positions
plot_df <- full_df %>%
  filter(condition %in% c("mel_only","time_5min","time_10min","time_20min","both_axes_10")) %>%
  mutate(
    subgroup = recode(condition,
                      "time_5min"    = "5 min",
                      "time_10min"   = "10 min",
                      "time_20min"   = "20 min",
                      "mel_only"     = "Mel Only",
                      "both_axes_10" = "Mel + Time"
    ),
    xpos = case_when(   # tighter spacing inside groups
      subgroup == "5 min"      ~ 1,
      subgroup == "10 min"     ~ 2,
      subgroup == "20 min"     ~ 3,
      subgroup == "Mel Only"   ~ 5,   # <- jump to create gap
      subgroup == "Mel + Time" ~ 6
    )
  )

# custom labels for x-axis
x_labels <- c(
  "5 min","10 min","20 min","Mel Only","Mel + Time"
)
#x_breaks <- c(1,2,3,5,6)
x_breaks <- c(1,1.5,2,3.5,4)

# assign distinct colors
cols <- c(
  "5 min"      = "#1b9e77",
  "10 min"     = "#d95f02",
  "20 min"     = "#7570b3",
  "Mel Only"   = "#e7298a",
  "Mel + Time" = "#66a61e"
)

ggplot(plot_df, aes(x = xpos, y = delta, fill = subgroup)) +
  geom_half_violin(side="l", trim=FALSE, alpha=0.6, width=0.6) +
  geom_boxplot(width=0.10, outlier.shape=NA,
               position=position_nudge(x=0.1), alpha=0.8) +
  geom_half_point(side="r", range_scale=0.4,
                  position=position_nudge(x=0.1),
                  alpha=0.4, size=1) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = c(0.05,0.05)) +
  scale_fill_manual(values=cols) +
  labs(
    title = "Raincloud Plot of ΔDLMO by Noise Condition",
    x = NULL,
    y = "Δ DLMO (hours)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x        = element_text(angle = 0, hjust = 1),
    legend.position    = "none",
    panel.grid.major.x = element_blank()
  ) +
  annotate("text", x = 2,   y = 2, label = "Time Only", fontface = "bold") +
  annotate("text", x = 5.5, y = 2, label = "Mel Group", fontface = "bold")


#######

library(dplyr)
library(ggplot2)
library(gghalves)

plot_df <- full_df %>%
  filter(condition %in% c("mel_only","time_5min","time_10min","time_20min","both_axes_10")) %>%
  mutate(
    subgroup = recode(condition,
                      "time_5min"    = "5 min",
                      "time_10min"   = "10 min",
                      "time_20min"   = "20 min",
                      "mel_only"     = "Mel Only",
                      "both_axes_10" = "Mel + Time"
    )
  )

# define order with a dummy "gap"
subgroup_levels <- c("5 min","10 min","20 min","", "Mel Only","Mel + Time")
plot_df <- plot_df %>% mutate(subgroup = factor(subgroup, levels = subgroup_levels))

# color palette
cols <- c(
  "5 min"      = "#1b9e77",
  "10 min"     = "#d95f02",
  "20 min"     = "#7570b3",
  "Mel Only"   = "#e7298a",
  "Mel + Time" = "#66a61e"
)

# nudge amount for boxplots/means
point_nudge <- 0.15

ggplot(plot_df, aes(x = subgroup, y = delta, fill = subgroup)) +
  # half violin on left
  geom_half_violin(side="l", alpha=0.6, width=0.9, scale="width", trim=FALSE) +
  # half boxplot on right
  geom_half_boxplot(side="r", outlier.shape=NA, width=0.2, color="black", fill=NA) +
  # individual points
  geom_point(
    aes(x = as.numeric(subgroup) + point_nudge, y = delta),
    shape = 21, size = 1.5, stroke = 0.2, color = "black", alpha = 0.6,
    position = position_jitter(width = 0.05, height = 0)
  ) +
  # mean marker
  stat_summary(
    fun = mean, geom = "point",
    shape = 21, size = 2.5, fill = "white", color = "black",
    position = position_nudge(x = point_nudge)
  ) +
  scale_fill_manual(values = cols, na.value = NA) +
  scale_x_discrete(
    labels = function(x) ifelse(x == "", "", x)  # hide dummy gap label
  ) +
  labs(
    title = "DLMO estimate sensitivty to noise",
    x = NULL,
    y = "Δ DLMO (hours)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position    = "none",
    panel.grid.major.x = element_blank()
  ) +
  annotate("text", x = 2,   y = 2, label = "Time Only", fontface = "bold") +
  annotate("text", x = 5.5, y = 2, label = "Mel Group", fontface = "bold")

######
library(dplyr)
library(ggplot2)
library(gghalves)

plot_df <- full_df %>%
  filter(condition %in% c("mel_only","time_5min","time_10min","time_20min","both_axes_10")) %>%
  mutate(
    subgroup = recode(condition,
                      "time_5min"    = "5 min",
                      "time_10min"   = "10 min",
                      "time_20min"   = "20 min",
                      "mel_only"     = "Mel Only",
                      "both_axes_10" = "Mel + 10 min"
    )
  )

# define order with dummy "gap"
subgroup_levels <- c("5 min","10 min","20 min","", "Mel Only","Mel + 10 min")
plot_df <- plot_df %>% mutate(subgroup = factor(subgroup, levels = subgroup_levels))

# color palette
cols <- c(
  "5 min"      = "#1b9e77",
  "10 min"     = "#d95f02",
  "20 min"     = "#7570b3",
  "Mel Only"   = "#e7298a",
  "Mel + 10 min" = "#66a61e"
)

# nudge amount
nudge_amt <- 0.2

ggplot(plot_df, aes(x = subgroup, y = delta, fill = subgroup)) +
  # half violin, shifted right
  geom_half_violin(side="l", alpha=0.6, width=0.9, scale="width", trim=FALSE,
                   position = position_nudge(x = -0.05)) +
  # half boxplot, also shifted right
  geom_half_boxplot(side="r", outlier.shape=NA, width=0.2, color="black", fill=NA,
                    position = position_nudge(x = 0.1)) +
  # dots stay centered (no nudge!)
  geom_point(
    shape = 21, size = 1.5, stroke = 0.2, color = "black", alpha = 0.6,
    position = position_jitter(width = 0.05, height = 0)
  ) +
  # mean marker, shifted right
  stat_summary(
    fun = mean, geom = "point",
    shape = 21, size = 2.5, fill = "white", color = "black",
    position = position_nudge(x = 0.1)
  ) +
  scale_fill_manual(values = cols, na.value = NA) +
  scale_x_discrete(labels = function(x) ifelse(x == "", "", x)) +
  labs(
    title = "DLMO estimate sensitivity to noise",
    x = NULL,
    y = "Δ DLMO (hours)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position    = "none",
    panel.grid.major.x = element_blank()
  ) +
  annotate("text", x = 2,   y = 2, label = "Time Only", fontface = "bold") +
  annotate("text", x = 4.5, y = 2, label = "Melatonin", fontface = "bold")

library(dplyr)

summary_stats <- plot_df %>%
  filter(subgroup != "") %>%  # drop the dummy gap
  group_by(subgroup) %>%
  summarise(
    mean_delta = mean(delta, na.rm = TRUE),
    sd_delta   = sd(delta, na.rm = TRUE),
    n          = n()
  )

print(summary_stats)
write.csv(summary_stats, "blume_dlmo_noise_summary_stats.csv", row.names = FALSE)


#####
# mean and SD plot
#####
# Line plot of mean ± SD
p_summary <- ggplot(summary_stats, aes(x = subgroup, y = mean_delta, group = 1)) +
  # Ribbon for ± SD
  geom_ribbon(
    aes(ymin = mean_delta - sd_delta, ymax = mean_delta + sd_delta),
    fill = "skyblue", alpha = 0.4
  ) +
  # Mean line
  geom_line(color = "black", size = 1) +
  # Mean points
  geom_point(shape = 21, size = 3, fill = "white", color = "black") +
  # X-axis categories (discrete here, so no scale_x_continuous)
  labs(
    title = "Mean Δ DLMO ± SD by Condition",
    x = "Condition",
    y = "Δ DLMO (hours)"
  ) +
  theme_minimal(base_size = 13)

p_summary


library(dplyr)
library(ggplot2)

# Subset: time-only conditions
time_summary <- summary_stats %>%
  filter(subgroup %in% c("5 min", "10 min", "20 min"))

# Subset: melatonin conditions
mel_summary <- summary_stats %>%
  filter(subgroup %in% c("Mel Only", "Mel + 10 min"))

# Plot time-only group
p_time <- ggplot(time_summary, aes(x = subgroup, y = mean_delta, group = 1)) +
  geom_ribbon(aes(ymin = mean_delta - sd_delta, ymax = mean_delta + sd_delta),
              fill = "skyblue", alpha = 0.4) +
  geom_line(color = "black", size = 1) +
  geom_point(shape = 21, size = 3, fill = "white", color = "black") +
  labs(
    title = "Time Only Conditions",
    x = "Condition",
    y = "Δ DLMO (hours)"
  ) +
  theme_minimal(base_size = 13)

# Plot melatonin group
p_mel <- ggplot(mel_summary, aes(x = subgroup, y = mean_delta, group = 1)) +
  geom_ribbon(aes(ymin = mean_delta - sd_delta, ymax = mean_delta + sd_delta),
              fill = "lightgreen", alpha = 0.4) +
  geom_line(color = "black", size = 1) +
  geom_point(shape = 21, size = 3, fill = "white", color = "black") +
  labs(
    title = "Melatonin Conditions",
    x = "Condition",
    y = "Δ DLMO (hours)"
  ) +
  theme_minimal(base_size = 13)

p_time
p_mel

library(patchwork)

(p_time | p_mel) +
  plot_annotation(
    title = "Mean Δ DLMO ± SD by Condition Groups",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

y_limits <- c(-1, 1)  # adjust as needed

p_time <- p_time + coord_cartesian(ylim = y_limits)
p_mel  <- p_mel  + coord_cartesian(ylim = y_limits)


p_mel <- p_mel +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  )

(p_time | p_mel) +
  plot_annotation(
    title = "Mean Δ DLMO ± SD by Condition Groups",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

