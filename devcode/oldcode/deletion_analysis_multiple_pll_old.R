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
library(lubridate)
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

decimal_to_posixct <- function(decimal_hour, reference_times) {
  start_time <- floor_date(min(reference_times), unit = "day")
  start_time + seconds(decimal_hour * 3600)
}

# -----------------------------
# multiple deletion analysis
# -----------------------------
run_multiple_deletion <- function(profile_id, df) {
  full_dlmo_result <- tryCatch({
    calculate_dlmo(df, threshold = 5)
  }, error = function(e) {
    message("Full DLMO failed for ", profile_id, ": ", e$message)
    return(NULL)
  })

  if (is.null(full_dlmo_result)) return(tibble())

  full_dlmo <- extract_dlmo_value(full_dlmo_result)
  full_dlmo_time <- decimal_to_posixct(full_dlmo, df$datetime)

  # percentages <- c(10, 20, 30, 40, 50) # change back to this!!!
  percentages <- c(10, 20)

  # Outer map per percentage, inner map per replicate
  multi_deletion <- map_dfr(percentages, function(pct) {
    n_del <- floor(pct / 100 * nrow(df))
    if (n_del < 1) return(NULL)

    future_map_dfr(1:2, function(rep) { #change this back to 100!!!
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
        deleted_timepoints = list(df$datetime[idx]),
        delta_dlmo = dlmo_deleted - full_dlmo,
        full_dlmo_time = full_dlmo_time,
        deleted_minutes_from_dlmo = list(relative_minutes_to_dlmo(df$datetime[idx], full_dlmo_time))
      )
    }, .options = furrr_options(seed = TRUE))
  })

  return(multi_deletion)
}

# -----------------------------
# Load and process all profiles
# -----------------------------
completed_ids <- list.files("multiple_deletion_results_blume", pattern = "\\.rds$") %>%
  file_path_sans_ext()

#profile_folder <- system.file("Documents/dlmoR/inst/extdata/", package = "dlmoR")
profile_folder <- "inst/extdata"
csv_files <- list.files(profile_folder, pattern = "\\.csv$", full.names = TRUE)

profiles <- set_names(csv_files, file_path_sans_ext(basename(csv_files))) %>%
  map(read_csv, show_col_types = FALSE)

profiles <- profiles[!names(profiles) %in% completed_ids]
profiles <- profiles[1:2]

# Enable nested parallelism: outer = multisession, inner = sequential
Sys.setenv("R_FUTURE_FORK_ENABLE" = "false")
options(future.rng.onMisuse = "ignore")
plan(list(multisession, sequential))
handlers(global = TRUE)

safe_run <- safely(run_multiple_deletion)

with_progress({
  p <- progressor(along = profiles)
  dir.create("deletion_partial_results_blume", showWarnings = FALSE)
  results_list <- future_imap(profiles, function(df, id) {
    result <- safe_run(id, df)
    if (!is.null(result$result)) {
      saveRDS(result$result, file = paste0("multiple_deletion_results_blume/", id, ".rds"))
    }
    p()
    result
  })
})

# -----------------------------
# Summary of processing
# -----------------------------
saved_files <- list.files("multiple_deletion_results_blume", pattern = "\\.rds$", full.names = TRUE)
all_results <- map_dfr(saved_files, readRDS)

num_processed <- length(results_list)
num_success <- sum(map_lgl(results_list, ~ !is.null(.x$result)))
num_failed <- num_processed - num_success

message("Summary:")
message("Profiles attempted: ", num_processed)
message("  - Successful: ", num_success)
message("  - Failed: ", num_failed)
message("Results saved to 'multiple_deletion_results_blume/'")

error_log <- imap_dfr(results_list, function(res, id) {
  if (!is.null(res$error)) {
    tibble(profile = id, error_message = res$error$message)
  }
})

write_csv(error_log, "multiple_deletion_errors_blume.csv")


# -----------------------------
# Deletion Raster Plot (by replicate, color by profile)
# -----------------------------
plot_deletion_raster_by_replicate <- function(results_df) {
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(forcats)

  deletions_long <- results_df %>%
    select(profile, percentage_deleted, replicate, deleted_minutes_from_dlmo) %>%
    unnest(deleted_minutes_from_dlmo) %>%
    rename(relative_minute = deleted_minutes_from_dlmo)

  deletions_long <- deletions_long %>%
    mutate(
      profile = as.factor(profile),
      replicate_id = interaction(profile, replicate, sep = ":"),
      replicate_id = fct_inorder(replicate_id)
    )

  ggplot(deletions_long, aes(x = relative_minute, y = replicate_id, fill = profile)) +
    geom_tile(height = 0.8) +
    facet_wrap(~percentage_deleted, scales = "free_y", ncol = 1) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    scale_x_continuous(
      breaks = seq(-800, 800, by = 60),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      title = "Deleted Timepoints by Replicate (Relative to DLMO)",
      subtitle = "Each tile = one deleted sample in one replicate",
      x = "Time relative to DLMO (minutes)",
      y = "Replicate",
      fill = "Profile"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "right"
    )
}
plot_deletion_raster_by_replicate(all_results)

# -----------------------------
# Deletion Violin Plots
# -----------------------------
plot_dlmo_deletion_violin <- function(results_df) {
  library(ggplot2)
  library(dplyr)
  library(gghalves)  # required for geom_half_violin, etc.

  results_df %>%
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

    # Half boxplot on the right
    geom_half_boxplot(
      side = "r",
      outlier.shape = NA,
      width = 0.2,
      color = "black",
      fill = NA
    ) +

    # Optional: Add mean point
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
}
plot_dlmo_deletion_violin(all_results)


