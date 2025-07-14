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
# Scenario 1 analysis only
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

  results <- map_dfr(seq_len(nrow(df)), function(i) {
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
  })

  return(results)
}

# -----------------------------
# Load and process profiles
# -----------------------------
profile_folder <- system.file("extdata/", package = "dlmoR")
csv_files <- list.files(profile_folder, pattern = "\\.csv$", full.names = TRUE)

profiles <- set_names(csv_files, file_path_sans_ext(basename(csv_files))) %>%
  map(read_csv, show_col_types = FALSE)

# Optional: test subset
# profiles <- profiles[1:2]

plan(multicore, workers = parallel::detectCores() - 1)
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
#
# plot1 <- all_results %>%
#   ggplot(aes(x = delta_minutes_from_dlmo, y = factor(profile), fill = delta_dlmo)) +
#   geom_tile(color = "white", linewidth = 0.2) +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Δ DLMO (h)") +
#   labs(
#     title = "Sensitivity of DLMO to Single-Point Deletions",
#     x = "Time of Deleted Sample (minutes relative to full-profile DLMO)",
#     y = "Profile"
#   ) +
#   theme_minimal()
#
# print(plot1)

# install.packages("forcats")
#
# library(forcats)
# library(scales)
# install.packages("RColorBrewer")
# library(RColorBrewer)
#
# # Compute symmetric color limits based on your data
# max_abs <- max(abs(all_results$delta_dlmo), na.rm = TRUE)
#
# dlmo_range <- range(all_results$delta_dlmo, na.rm = TRUE)
# dlmo_range <- c(floor(dlmo_range[1]), ceiling(dlmo_range[2]))
#
# # Compute global range for limits
# dlmo_range <- range(all_results$delta_dlmo, na.rm = TRUE)
#
# plot1 <- all_results %>%
#   mutate(
#     profile = forcats::fct_rev(factor(profile)),
#     rounded_minutes = round(delta_minutes_from_dlmo)
#   ) %>%
#   ggplot(aes(x = rounded_minutes, y = profile, fill = delta_dlmo)) +
#   geom_tile(color = NA, height = 0.5) +
#   geom_tile(color = NA, height = 0.9, width = 60)+
#   geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5)+
#   scale_fill_gradient2(
#     low = "#4575b4", mid = "#f7f7f7", high = "#d73027",
#     midpoint = 0, limits = dlmo_range,
#     name = "Δ DLMO (h)"
#   ) +
#   scale_x_continuous(
#     breaks = seq(-800, 800, by = 120),  # adjust as needed based on your time range
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
# # Get all profile names
# profile_names <- unique(all_results$profile)
#
# # Exclude the last one
# profiles_to_plot <- head(profile_names, -1)
#
# # Filter data for plotting
# filtered_results <- all_results %>%
#   filter(profile %in% profiles_to_plot)
#
#
# # Compute global range for limits
# # Compute global range for limits
# dlmo_range <- range(all_results$delta_dlmo, na.rm = TRUE)
#
# # Get all profile names and remove the last one
# profile_names <- unique(all_results$profile)
# profiles_to_plot <- head(profile_names, -1)
#
# # Filter data for plotting
# filtered_results <- all_results %>%
#   filter(profile %in% profiles_to_plot) %>%
#   mutate(
#     profile = forcats::fct_rev(factor(profile)),
#     rounded_minutes = round(delta_minutes_from_dlmo)
#   )
#
# # Plot
# plot1 <- ggplot(data = filtered_results, aes(x = rounded_minutes, y = profile, fill = delta_dlmo)) +
#   geom_tile(color = NA, height = 0.5) +
#   geom_tile(color = NA, height = 0.9, width = 60) +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
#   scale_fill_gradient2(
#     low = "#4575b4", mid = "#f7f7f7", high = "#d73027",
#     midpoint = 0, limits = dlmo_range,
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

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(forcats)
library(scico)

# Get all profile names and remove the last one
profile_names <- unique(all_results$profile)
profiles_to_plot <- head(profile_names, -1)

# Compute global color scale limits (symmetric)
dlmo_range <- range(all_results$delta_dlmo, na.rm = TRUE)
dlmo_range <- c(-max(abs(dlmo_range)), max(abs(dlmo_range)))

# Filter data for plotting
filtered_results <- all_results %>%
  filter(profile %in% profiles_to_plot) %>%
  mutate(
    profile = forcats::fct_rev(factor(profile)),
    rounded_minutes = round(delta_minutes_from_dlmo)
  )

# Plot
plot1 <- ggplot(data = filtered_results, aes(x = rounded_minutes, y = profile, fill = delta_dlmo)) +
  geom_tile(color = NA, height = 0.5) +
  geom_tile(color = NA, height = 0.9, width = 60) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  scale_fill_scico(
    palette = "vik",       # perceptually uniform diverging palette
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

# Print the plot
print(plot1)

