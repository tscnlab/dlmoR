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
# Scenario 2 analysis only
# -----------------------------
run_scenario2_only <- function(profile_id, df) {
  full_dlmo_result <- tryCatch({
    calculate_dlmo(df, threshold = 5)
  }, error = function(e) {
    message("Full DLMO failed for ", profile_id, ": ", e$message)
    return(NULL)
  })

  if (is.null(full_dlmo_result)) return(tibble())

  full_dlmo <- extract_dlmo_value(full_dlmo_result)
  full_dlmo_time <- decimal_to_posixct(full_dlmo, df$datetime)

  percentages <- c(10, 20, 30, 40, 50)

  # Outer map per percentage, inner map per replicate
  scenario2 <- map_dfr(percentages, function(pct) {
    n_del <- floor(pct / 100 * nrow(df))
    if (n_del < 1) return(NULL)

    future_map_dfr(1:100, function(rep) {
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

  return(scenario2)
}

# -----------------------------
# Load and process all profiles
# -----------------------------
completed_ids <- list.files("deletion_partial_results_blume", pattern = "\\.rds$") %>%
  file_path_sans_ext()

#profile_folder <- system.file("Documents/dlmoR/inst/extdata/", package = "dlmoR")
profile_folder <- "inst/extdata"
csv_files <- list.files(profile_folder, pattern = "\\.csv$", full.names = TRUE)

profiles <- set_names(csv_files, file_path_sans_ext(basename(csv_files))) %>%
  map(read_csv, show_col_types = FALSE)

profiles <- profiles[!names(profiles) %in% completed_ids]
profiles <- profiles[1:40]

# Enable nested parallelism: outer = multisession, inner = sequential
Sys.setenv("R_FUTURE_FORK_ENABLE" = "false")
options(future.rng.onMisuse = "ignore")
plan(list(multisession, sequential))
handlers(global = TRUE)

safe_run <- safely(run_scenario2_only)

with_progress({
  p <- progressor(along = profiles)
  dir.create("deletion_partial_results_blume", showWarnings = FALSE)
  results_list <- future_imap(profiles, function(df, id) {
    result <- safe_run(id, df)
    if (!is.null(result$result)) {
      saveRDS(result$result, file = paste0("deletion_partial_results_blume/", id, ".rds"))
    }
    p()
    result
  })
})

# -----------------------------
# Summary of processing
# -----------------------------
saved_files <- list.files("deletion_partial_results_blume", pattern = "\\.rds$", full.names = TRUE)
all_results <- map_dfr(saved_files, readRDS)

num_processed <- length(results_list)
num_success <- sum(map_lgl(results_list, ~ !is.null(.x$result)))
num_failed <- num_processed - num_success

message("Summary:")
message("Profiles attempted: ", num_processed)
message("  - Successful: ", num_success)
message("  - Failed: ", num_failed)
message("Partial results saved to 'deletion_partial_results_blume/'")

error_log <- imap_dfr(results_list, function(res, id) {
  if (!is.null(res$error)) {
    tibble(profile = id, error_message = res$error$message)
  }
})

write_csv(error_log, "dlmo_deletion_errors_blume.csv")


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
