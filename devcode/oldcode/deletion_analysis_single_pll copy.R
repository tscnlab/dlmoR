# # -----------------------------
# # Required packages
# # -----------------------------
# library(readr)
# library(dplyr)
# library(purrr)
# library(tibble)
# library(tools)
# library(progressr)
# library(furrr)
# library(future)
# library(lubridate)
# library(ggplot2)
# library(forcats)
# library(scico)
#
# # -----------------------------
# # Utility functions
# # -----------------------------
# extract_dlmo_value <- function(dlmo_result) {
#   dlmo_result$ip$inflection_point_fine$x
# }
#
# relative_minutes_to_dlmo <- function(timestamps, dlmo_time) {
#   as.numeric(difftime(timestamps, dlmo_time, units = "mins"))
# }
#
# decimal_to_posixct <- function(decimal_hour, reference_times) {
#   start_time <- floor_date(min(reference_times), unit = "day")
#   start_time + seconds(decimal_hour * 3600)
# }
#
# # -----------------------------
# # Scenario 1 analysis only (parallelized internally)
# # -----------------------------
# run_scenario1_only <- function(profile_id, df) {
#   full_dlmo_result <- tryCatch({
#     calculate_dlmo(df, threshold = 5)
#   }, error = function(e) {
#     message("Full DLMO failed for ", profile_id, ": ", e$message)
#     return(NULL)
#   })
#
#   if (is.null(full_dlmo_result)) return(tibble())
#
#   full_dlmo <- extract_dlmo_value(full_dlmo_result)
#   full_dlmo_time <- decimal_to_posixct(full_dlmo_result$ip$inflection_point_fine$x, full_dlmo_result$prof$datetime)
#
#   results <- future_map_dfr(seq_len(nrow(df)), function(i) {
#     df_deleted <- df[-i, ]
#     dlmo_deleted <- tryCatch({
#       extract_dlmo_value(calculate_dlmo(df_deleted, threshold = 5))
#     }, error = function(e) NA_real_)
#
#     tibble(
#       profile = profile_id,
#       deleted_time = df$datetime[i],
#       delta_minutes_from_dlmo = relative_minutes_to_dlmo(df$datetime[i], full_dlmo_time),
#       delta_dlmo = dlmo_deleted - full_dlmo
#     )
#   }, .options = furrr_options(seed = TRUE))
#
#   return(results)
# }
#
# # -----------------------------
# # Load and process profiles
# # -----------------------------
# profile_folder <- system.file("extdata/", package = "dlmoR")
# csv_files <- list.files(profile_folder, pattern = "\\.csv$", full.names = TRUE)
#
# # Identify processed profiles
# result_files <- list.files("deletion_scenario1_results", pattern = "\\.rds$")
# processed_ids <- file_path_sans_ext(basename(result_files))
#
# # Load only unprocessed profiles
# profiles <- set_names(csv_files, file_path_sans_ext(basename(csv_files)))
# profiles <- profiles[!names(profiles) %in% processed_ids]
# profiles <- map(profiles, read_csv, show_col_types = FALSE)
#
# # Parallel plan: nested parallelization
# plan(nesting(multisession, multisession))
# handlers(global = TRUE)
#
# safe_run <- safely(run_scenario1_only)
#
# with_progress({
#   p <- progressor(along = profiles)
#   dir.create("deletion_scenario1_results", showWarnings = FALSE)
#   results_list <- future_imap(profiles, function(df, id) {
#     result <- safe_run(id, df)
#     if (!is.null(result$result)) {
#       saveRDS(result$result, file = file.path("deletion_scenario1_results", paste0(id, ".rds")))
#     }
#     p()
#     result
#   })
# })
#
# # -----------------------------
# # Combine and plot results
# # -----------------------------
# saved_files <- list.files("deletion_scenario1_results", pattern = "\\.rds$", full.names = TRUE)
# all_results <- map_dfr(saved_files, readRDS)
#
# # Prepare data for plotting
# profile_names <- unique(all_results$profile)
# profiles_to_plot <- head(profile_names, -1)
#
# dlmo_range <- range(all_results$delta_dlmo, na.rm = TRUE)
# dlmo_range <- c(-max(abs(dlmo_range)), max(abs(dlmo_range)))
#
# filtered_results <- all_results %>%
#   filter(profile %in% profiles_to_plot) %>%
#   mutate(
#     profile = forcats::fct_rev(factor(profile)),
#     rounded_minutes = round(delta_minutes_from_dlmo)
#   )
#
# plot1 <- ggplot(data = filtered_results, aes(x = rounded_minutes, y = profile, fill = delta_dlmo)) +
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
# energy
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
run_scenario1_only <- function(profile_id, df) {
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
result_files <- list.files("deletion_scenario1_results", pattern = "\\.rds$")
processed_ids <- file_path_sans_ext(basename(result_files))

# Load only unprocessed profiles
profiles <- set_names(csv_files, file_path_sans_ext(basename(csv_files)))
profiles <- profiles[!names(profiles) %in% processed_ids]
profiles <- map(profiles, read_csv, show_col_types = FALSE)

# Parallel plan: single-level multisession for M4 chip
Sys.setenv("R_FUTURE_FORK_ENABLE" = "false")
options(future.rng.onMisuse = "ignore")
plan(multisession, workers = parallel::detectCores(logical = FALSE))
handlers(global = TRUE)

safe_run <- safely(run_scenario1_only)

with_progress({
  p <- progressor(along = profiles)
  dir.create("deletion_scenario1_results", showWarnings = FALSE)
  results_list <- future_imap(profiles, function(df, id) {
    result <- safe_run(id, df)
    if (!is.null(result$result)) {
      saveRDS(result$result, file = file.path("deletion_scenario1_results", paste0(id, ".rds")))
    }
    p()
    result
  })
})

# -----------------------------
# Combine and plot results
# -----------------------------
saved_files <- list.files("deletion_scenario1_results", pattern = "\\.rds$", full.names = TRUE)
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
    na.value = "grey40",  # or "lightgrey", "black", "red", etc.
    name = "Δ DLMO (h)"
  )+
  scale_x_continuous(
    breaks = seq(-800, 800, by = 60),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Sensitivity of DLMO to single-point deletions",
    subtitle = "Blume et al. dataset",
    x = "Deleted sample time (minutes relative to full-profile DLMO)",
    y = "Melatonin profile"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    panel.grid = element_blank()
  )+
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

print(plot1)

# line plot
library(scales)  # for rescale

# Prepare summary
summary_df <- filtered_results %>%
  group_by(rounded_minutes) %>%
  summarise(
    n = sum(!is.na(delta_dlmo)),
    mad_dlmo = mad(delta_dlmo, na.rm = TRUE)
  )

# Rescale n to align with y-axis for overlay
summary_df <- summary_df %>%
  mutate(n_scaled = rescale(n, to = c(0, max(mad_dlmo, na.rm = TRUE))))

# Plot
library(dplyr)
library(ggplot2)

plot1a <- filtered_results %>%
  mutate(rounded_minutes = round(delta_minutes_from_dlmo / 10) * 10) %>%  # <-- change bin size here
  group_by(rounded_minutes) %>%
  summarise(
    n = sum(!is.na(delta_dlmo)),
    mad_dlmo = mad(delta_dlmo, na.rm = TRUE)
  ) %>%
  ggplot(aes(x = rounded_minutes, y = mad_dlmo)) +
  geom_line(color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    title = "DLMO sensitivity to deletion at each timepoint",
    x = "Deleted Sample Time (min relative to DLMO)",
    y = "MAD of Δ DLMO (hours)"
  ) +
  theme_minimal()

plot(plot1a)

#violin plot
plot1b <- filtered_results %>%
  mutate(time_bin = round(delta_minutes_from_dlmo / 10) * 10) %>%
  ggplot(aes(x = factor(time_bin), y = delta_dlmo)) +
  geom_violin(
    aes(fill = factor(time_bin)),
    scale = "width",
    adjust = 0.7,
    trim = TRUE,
    color = NA,
    alpha = 0.6
  ) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 1.5, fill = "white") +
  labs(
    title = "Sensitivity of DLMO to deletion time (across all profiles)",
    x = "Deleted Time (min relative to DLMO, binned)",
    y = "Δ DLMO (hours)"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")+
  scale_x_discrete(
  breaks = function(x) {
    x_vals <- as.numeric(as.character(x))
    breaks <- x_vals[which(x_vals %% 20 == 0)]
    unique(as.character(breaks))
  }
)

plot(plot1b)

# raincloud plot
library(dplyr)
library(ggplot2)
library(gghalves)

# Prepare data
binned_data <- filtered_results %>%
  filter(!is.na(delta_dlmo), !is.na(rounded_minutes)) %>%
  mutate(
    binned_minutes = floor(rounded_minutes / 10) * 10,
    binned_minutes = factor(binned_minutes)  # ensure discrete x
  )

# Filter bins with at least 2 values (geom_half_violin needs this)
binned_counts <- binned_data %>%
  group_by(binned_minutes) %>%
  filter(n() >= 2) %>%
  ungroup()

# X-axis breaks: every 20 min
breaks_to_show <- levels(binned_counts$binned_minutes)[as.numeric(as.character(levels(binned_counts$binned_minutes))) %% 20 == 0]

# Plot
plot1raincloud <- ggplot(binned_counts, aes(x = binned_minutes, y = delta_dlmo)) +
  geom_half_violin(
    side = "l",
    fill = "skyblue",
    alpha = 0.6,
    width = 1,
    scale = "width",
    trim = TRUE,
    bw = 0.2,
    adjust = 0.5,
    color = NA
  ) +
  geom_half_point(
    side = "r",
    shape = 21,
    fill = "black",
    color = "black",
    alpha = 0.4,
    size = 1.2,
    width = 0.2
  ) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 1.8, fill = "white") +
  scale_x_discrete(breaks = breaks_to_show) +
  labs(
    title = "Raincloud Plot: Δ DLMO by Deletion Time (binned)",
    x = "Deleted Time (min relative to DLMO)",
    y = "Δ DLMO (h)"
  ) +
  theme_minimal(base_size = 13)

print(plot1raincloud)

plot1raincloud <- ggplot(binned_counts, aes(x = binned_minutes, y = delta_dlmo)) +

  # Half violin (distribution)
  geom_half_violin(
    side = "l",
    fill = "skyblue",
    alpha = 0.6,
    width = 1,
    scale = "width",
    trim = TRUE,
    bw = 0.2,
    adjust = 0.5,
    color = NA
  ) +

  # Box plot in the middle
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    alpha = 0.6
  ) +

  # Jittered points (rain)
  geom_half_point(
    side = "r",
    shape = 21,
    fill = "black",
    color = "black",
    alpha = 0.2,  # << lighter for readability
    size = 1,
    width = 0.2
  ) +

  scale_x_discrete(breaks = breaks_to_show) +

  labs(
    title = "DLMO sensitivity to deletion time",
    subtitle = "Blume et al. dataset",
    x = "Deleted timepoint (minutes relative to full-profile DLMO)",
    y = "Δ DLMO (h)"
  ) +
  theme_minimal(base_size = 13)
print(plot1raincloud)

# ribbon plot
library(dplyr)
library(ggplot2)

# Bin and summarize
summary_data <- filtered_results %>%
  mutate(binned_minutes = floor(rounded_minutes / 10) * 10) %>%
  group_by(binned_minutes) %>%
  summarise(
    mean_dlmo = mean(delta_dlmo, na.rm = TRUE),
    sd_dlmo = sd(delta_dlmo, na.rm = TRUE),
    n = sum(!is.na(delta_dlmo)),
    sem_dlmo = sd_dlmo / sqrt(n)
  )

# Plot: mean ± SD ribbon (or change to SEM)
plot1ribbon <- ggplot(summary_data, aes(x = binned_minutes, y = mean_dlmo)) +
  geom_ribbon(aes(ymin = mean_dlmo - sd_dlmo, ymax = mean_dlmo + sd_dlmo), fill = "skyblue", alpha = 0.4) +
  geom_line(color = "black", size = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    title = "Systematic DLMO Bias from Deleting Single Timepoints",
    x = "Deleted Sample Time (min relative to DLMO)",
    y = "Mean Δ DLMO (hours)"
  ) +
  theme_minimal(base_size = 13)
plot(plot1ribbon)


library(dplyr)
library(ggplot2)

# Dummy variable to enable legend mapping
summary_data$group <- "Bias (Mean ± SD)"

plot1ribbon <- ggplot(summary_data, aes(x = binned_minutes)) +

  # Ribbon with manual fill legend
  geom_ribbon(
    aes(ymin = mean_dlmo - sd_dlmo, ymax = mean_dlmo + sd_dlmo, fill = group),
    alpha = 0.4
  ) +

  # Smoothed mean line with manual color legend
  geom_smooth(
    aes(y = mean_dlmo, color = group),
    method = "loess", span = 0.25, se = FALSE,
    linewidth = 1.2
  ) +

  # Reference lines
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray60") +

  # Legend and color customization
  scale_fill_manual(name = NULL, values = c("Bias (Mean ± SD)" = "skyblue")) +
  scale_color_manual(name = NULL, values = c("Bias (Mean ± SD)" = "black")) +

  labs(
    title = "Systematic DLMO Bias from Deleting Single Timepoints",
    x = "Deleted Sample Time (min relative to DLMO)",
    y = "Mean Δ DLMO (hours)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.justification = "left"
  )

plot(plot1ribbon)


plot_smoothed_bias <- ggplot(summary_data, aes(x = binned_minutes)) +
  geom_ribbon(
    aes(ymin = mean_dlmo - sd_dlmo, ymax = mean_dlmo + sd_dlmo, fill = group),
    alpha = 0.4
  ) +
  geom_smooth(
    aes(y = mean_dlmo, color = group),
    method = "loess", span = 0.15, se = FALSE,
    linewidth = 1
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray60") +
  scale_fill_manual(name = NULL, values = c("Bias (Mean ± SD)" = "skyblue")) +
  scale_color_manual(name = NULL, values = c("Bias (Mean ± SD)" = "black")) +
  labs(
    title = "DLMO sensitivity to deletion time",
    subtitle = "Blume et al. dataset",
    x = "Deleted timepoint (minutes relative to full-profile DLMO)",
    y = "Δ DLMO (h)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = c(0.95, 0.05),
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = "white", color = "gray80"),
    legend.title = element_blank(),
  )+
  theme(legend.background = element_blank())

plot(plot_smoothed_bias)

#
# plot1 <- ggplot(data = filtered_results, aes(x = rounded_minutes, y = profile)) +
#   # First: continuous delta_dlmo layer
#   geom_tile(
#     data = filter(filtered_results, !is.na(delta_dlmo)),
#     aes(fill = delta_dlmo),
#     color = NA, height = 0.9, width = 60
#   ) +
#   scale_fill_scico(
#     palette = "vik",
#     midpoint = 0,
#     limits = dlmo_range,
#     name = "Δ DLMO (h)"
#   ) +
#
#   # Switch fill scale
#   ggnewscale::new_scale_fill() +
#
#   # Now: NA values (discrete fill)
#   geom_tile(
#     data = filter(filtered_results, is.na(delta_dlmo)),
#     aes(fill = "Missing DLMO"),
#     color = NA, height = 0.9, width = 60
#   ) +
#   scale_fill_manual(
#     values = c("Missing DLMO" = "gray50"),
#     name = NULL,
#     guide = guide_legend(override.aes = list(color = NA))
#   ) +
#
#   # Axis, theme, etc.
#   geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
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

