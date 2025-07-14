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

# -----------------------------
# Utility functions
# -----------------------------
extract_dlmo_value <- function(dlmo_result) {
  dlmo_result$ip$inflection_point_fine$x
}

relative_minutes_to_dlmo <- function(timestamps, dlmo_time) {
  as.numeric(difftime(timestamps, dlmo_time, units = "mins"))
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

  percentages <- c(10, 20, 30, 40, 50)

  # Outer map per percentage, inner future_map for replicates
  scenario2 <- future_map_dfr(percentages, function(pct) {
    n_del <- floor(pct / 100 * nrow(df))
    if (n_del < 1) return(NULL)

    map_dfr(1:100, function(rep) {
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
        delta_dlmo = dlmo_deleted - full_dlmo
      )
    })
  }, .options = furrr_options(seed = TRUE))

  return(scenario2)
}

# -----------------------------
# Load and process all profiles
# -----------------------------
completed_ids <- list.files("deletion_partial_results_test2", pattern = "\\.rds$") %>%
  file_path_sans_ext()

profile_folder <- system.file("extdata/", package = "dlmoR")
csv_files <- list.files(profile_folder, pattern = "\\.csv$", full.names = TRUE)

profiles <- set_names(csv_files, file_path_sans_ext(basename(csv_files))) %>%
  map(read_csv, show_col_types = FALSE)

profiles <- profiles[!names(profiles) %in% completed_ids]
profiles <- profiles[1:2]

# Enable nested parallelism: outer = multisession, inner = sequential
Sys.setenv("R_FUTURE_FORK_ENABLE" = "false")
options(future.rng.onMisuse = "ignore")
plan(list(multisession, multisession))
handlers(global = TRUE)

safe_run <- safely(run_scenario2_only)

with_progress({
  p <- progressor(along = profiles)
  dir.create("deletion_partial_results_test2", showWarnings = FALSE)
  results_list <- future_imap(profiles, function(df, id) {
    result <- safe_run(id, df)
    if (!is.null(result$result)) {
      saveRDS(result$result, file = paste0("deletion_partial_results_test2/", id, ".rds"))
    }
    p()
    result
  })
})

# -----------------------------
# Summary of processing
# -----------------------------
saved_files <- list.files("deletion_partial_results_test2", pattern = "\\.rds$", full.names = TRUE)
all_results <- map_dfr(saved_files, readRDS)

num_processed <- length(results_list)
num_success <- sum(map_lgl(results_list, ~ !is.null(.x$result)))
num_failed <- num_processed - num_success

message("Summary:")
message("Profiles attempted: ", num_processed)
message("  - Successful: ", num_success)
message("  - Failed: ", num_failed)
message("Partial results saved to 'deletion_partial_results_test2/'")

error_log <- imap_dfr(results_list, function(res, id) {
  if (!is.null(res$error)) {
    tibble(profile = id, error_message = res$error$message)
  }
})

write_csv(error_log, "dlmo_deletion_errors_test2.csv")


plot2 <- all_results %>%
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
  geom_half_boxplot(
    side = "r",
    outlier.shape = NA,
    width = 0.2,
    color = "black",
    fill = NA
  )+
  # Optional: Add mean point or boxplot
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

print(plot2)


library(tidyr)

unnested_deleted <- all_results %>%
  select(profile, percentage_deleted, replicate, deleted_timepoints) %>%
  unnest(deleted_timepoints)

print(unnested_deleted, n = 180)
