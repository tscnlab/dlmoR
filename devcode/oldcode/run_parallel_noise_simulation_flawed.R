library(dplyr)
library(purrr)
library(dlmoR)
library(readr)
library(furrr)
library(future)
library(progressr)

# ────────────────────────────────────────────────────
# SETUP PARALLEL BACKEND
# ────────────────────────────────────────────────────
plan(multicore, workers = 64)
handlers(global = TRUE)

# ────────────────────────────────────────────────────
# LOAD PROFILES
# ────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
profile_folder <- if (length(args) > 0) args[[1]] else "inst/extdata"
csv_files <- list.files(profile_folder, pattern = "\\.csv$", full.names = TRUE)

out_dir <- "outputs/noisy_dlmo_results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

result_files  <- list.files(out_dir, pattern = "\\.rds$", full.names = TRUE)
processed_ids <- tools::file_path_sans_ext(basename(result_files))

profiles <- set_names(csv_files, tools::file_path_sans_ext(basename(csv_files)))
profiles_to_run <- profiles[!names(profiles) %in% processed_ids]

# ────────────────────────────────────────────────────
# NOISE PARAMETERS
# ────────────────────────────────────────────────────
cv_intra      <- 7.9
n_rep_default <- 20

enforce_min_gap <- function(times_dec, min_gap = 1/60) {
  for (i in seq(2, length(times_dec))) {
    if ((times_dec[i] - times_dec[i - 1]) < min_gap) {
      times_dec[i] <- times_dec[i - 1] + min_gap
    }
  }
  times_dec
}

# ────────────────────────────────────────────────────
# CORE FUNCTION (SINGLE PROFILE PROCESSING)
# ────────────────────────────────────────────────────
process_single_profile <- function(profile_id, path, n_rep = n_rep_default, cv = cv_intra) {
  df <- read_csv(path, show_col_types = FALSE) %>%
    mutate(datetime = as.POSIXct(datetime, tz = "UTC"))

  dec_t <- posixct_to_decimal(df$datetime, df$datetime[3])
  interp_fun <- splinefun(dec_t, df$melatonin, method = "monoH.FC")

  error_rows <- list()
  all_results <- list()
  all_clipping_logs <- list()

  run_one <- function(time_sd, add_mel_noise, condition, replicate) {
    t_pert <- enforce_min_gap(sort(dec_t + rnorm(length(dec_t), 0, time_sd)))
    mel_int_raw <- interp_fun(t_pert)

    clipping_log <- tibble()

    mel_int_neg_idx <- which(mel_int_raw < 0)
    if (length(mel_int_neg_idx) > 0) {
      clipping_log <- add_row(clipping_log, tibble(
        profile = profile_id,
        condition = condition,
        replicate = replicate,
        type = "interpolation",
        n_replaced = length(mel_int_neg_idx),
        values = paste(round(mel_int_raw[mel_int_neg_idx], 3), collapse = ";")
      ))
    }
    mel_int <- pmax(mel_int_raw, 0)

    if (add_mel_noise) {
      sd_vec <- (cv / 100) * mel_int
      mel_pert_raw <- mel_int + rnorm(length(mel_int), 0, sd_vec)

      mel_pert_neg_idx <- which(mel_pert_raw < 0)
      if (length(mel_pert_neg_idx) > 0) {
        clipping_log <- add_row(clipping_log, tibble(
          profile = profile_id,
          condition = condition,
          replicate = replicate,
          type = "perturbed",
          n_replaced = length(mel_pert_neg_idx),
          values = paste(round(mel_pert_raw[mel_pert_neg_idx], 3), collapse = ";")
        ))
      }

      mel_pert <- pmax(mel_pert_raw, 0)

      mean_mel_sd <- mean(sd_vec)
      min_mel_sd <- min(sd_vec)
      max_mel_sd <- max(sd_vec)
    } else {
      mel_pert <- mel_int
      mean_mel_sd <- min_mel_sd <- max_mel_sd <- 0
    }

    sim_df <- tibble(
      datetime  = decimal_to_posixct(t_pert, df$datetime[3]),
      melatonin = mel_pert
    )

    result <- tryCatch(calculate_dlmo(sim_df, threshold = 5), error = function(e) e)

    if (inherits(result, "error")) {
      error_rows <<- append(error_rows, list(tibble(
        profile = profile_id,
        condition = condition,
        replicate = replicate,
        error_message = result$message,
        sim_df = list(sim_df)
      )))
      return(list(result = NULL, log = clipping_log))
    }

    result_row <- tibble(
      profile = profile_id,
      sim_df = list(sim_df),
      time_sd = time_sd,
      mean_mel_sd = mean_mel_sd,
      min_mel_sd = min_mel_sd,
      max_mel_sd = max_mel_sd,
      dlmo_est = result$ip$inflection_point_fine$x,
      dlmo_full_result = list(result),
      condition = condition,
      replicate = replicate
    )

    list(result = result_row, log = clipping_log)
  }

  conditions <- list(
    clean = list(time_sd = 0, add_mel = FALSE),
    mel_only = list(time_sd = 0, add_mel = TRUE),
    time_5min = list(time_sd = 5 / 60, add_mel = FALSE),
    time_10min = list(time_sd = 10 / 60, add_mel = FALSE),
    time_20min = list(time_sd = 20 / 60, add_mel = FALSE),
    both_axes_10 = list(time_sd = 10 / 60, add_mel = TRUE)
  )

  for (cond in names(conditions)) {
    p <- conditions[[cond]]
    reps <- if (cond == "clean") 1 else n_rep

    for (i in seq_len(reps)) {
      run <- run_one(p$time_sd, p$add_mel, cond, i)
      if (!is.null(run$result)) {
        all_results <- append(all_results, list(run$result))
      }
      if (nrow(run$log) > 0) {
        all_clipping_logs <- append(all_clipping_logs, list(run$log))
      }
    }
  }

  if (length(all_clipping_logs) > 0) {
    write_csv(bind_rows(all_clipping_logs), file.path(out_dir, paste0(profile_id, "_clipping_log.csv")))
  }

  results_df <- bind_rows(all_results)

  clean_dlmo_est <- results_df %>%
    filter(condition == "clean", replicate == 1) %>%
    pull(dlmo_est)

  if (length(clean_dlmo_est) == 0 || is.na(clean_dlmo_est)) {
    warning("No valid clean DLMO estimate for profile: ", profile_id)
    clean_dlmo_est <- NA_real_
  }

  results_df <- results_df %>%
    mutate(
      clean_dlmo_ref = clean_dlmo_est,
      dlmo_error = clean_dlmo_ref - dlmo_est
    )

  list(
    results = results_df,
    errors = bind_rows(error_rows)
  )
}


# ────────────────────────────────────────────────────
# ERROR LOG CSV
# ────────────────────────────────────────────────────
error_log <- keep(results, ~ !is.null(.x$error)) %>%
  map_dfr(~ tibble(profile = .x$profile, error_message = .x$error$message))

write_csv(error_log, file.path(out_dir, "noise_simulation_errors.csv"))

# ────────────────────────────────────────────────────
# COMBINE ALL ITERATION-LEVEL ERRORS INTO ONE FILE
# ────────────────────────────────────────────────────
error_rds_files <- list.files(out_dir, pattern = "_iteration_errors\\.rds$", full.names = TRUE)

all_iteration_errors <- purrr::map_dfr(error_rds_files, readRDS)

readr::write_csv(all_iteration_errors, file.path(out_dir, "all_iteration_errors.csv"))

clipping_csvs <- list.files(out_dir, pattern = "_clipping_log\\.csv$", full.names = TRUE)

all_clipping_logs <- map_dfr(clipping_csvs, read_csv)

write_csv(all_clipping_logs, file.path(out_dir, "all_clipping_events.csv"))


# ────────────────────────────────────────────────────
# COMBINED RESULTS
# ────────────────────────────────────────────────────
# all_result_files <- list.files(out_dir, pattern = "\\.rds$", full.names = TRUE)
# all_results <- map_dfr(all_result_files, readRDS)

# saveRDS(all_results, file.path(out_dir, "dlmo_noise_all_results.rds"))

# # ────────────────────────────────────────────────────
# # NOISE PARAMETERS
# # ────────────────────────────────────────────────────
# cv_intra      <- 7.9
# n_rep_default <- 20
#
# enforce_min_gap <- function(times_dec, min_gap = 1/60) {
#   for (i in seq(2, length(times_dec))) {
#     if ((times_dec[i] - times_dec[i - 1]) < min_gap) {
#       times_dec[i] <- times_dec[i - 1] + min_gap
#     }
#   }
#   times_dec
# }
#
# # ────────────────────────────────────────────────────
# # CORE FUNCTION (SINGLE PROFILE PROCESSING)
# # ────────────────────────────────────────────────────
# process_single_profile <- function(profile_id, path, n_rep = n_rep_default, cv = cv_intra) {
#   df <- readr::read_csv(path, show_col_types = FALSE) %>%
#     dplyr::mutate(datetime = as.POSIXct(datetime, tz = "UTC"))
#
#   clean_res <- calculate_dlmo(df, threshold = 5)
#   dec_t <- posixct_to_decimal(df$datetime, df$datetime)
#   interp_fun <- splinefun(dec_t, df$melatonin, method = "monoH.FC")
#
#   time_sd_5  <- 5 / 60
#   time_sd_10 <- 10 / 60
#   time_sd_20 <- 20 / 60
#
#   error_rows <- list()
#
#   run_one <- function(time_sd, add_mel_noise, condition, replicate) {
#     t_pert <- enforce_min_gap(sort(dec_t + rnorm(length(dec_t), 0, time_sd)))
#     mel_int <- interp_fun(t_pert)
#
#     if (add_mel_noise) {
#       sd_vec <- (cv / 100) * mel_int
#       mel_pert <- mel_int + rnorm(length(mel_int), 0, sd_vec)
#       mean_mel_sd <- mean(sd_vec)
#       min_mel_sd <- min(sd_vec)
#       max_mel_sd <- max(sd_vec)
#     } else {
#       mel_pert <- mel_int
#       mean_mel_sd <- min_mel_sd <- max_mel_sd <- 0
#     }
#
#     sim_df <- tibble::tibble(
#       datetime  = decimal_to_posixct(t_pert, df$datetime),
#       melatonin = mel_pert
#     )
#
#     result <- tryCatch(calculate_dlmo(sim_df, threshold = 5), error = function(e) e)
#
#     if (inherits(result, "error")) {
#       error_rows <<- append(error_rows, list(tibble::tibble(
#         profile = profile_id,
#         condition = condition,
#         replicate = replicate,
#         error_message = result$message,
#         sim_df = list(sim_df)
#       )))
#       return(NULL)
#     }
#
#
#     tibble::tibble(
#       profile = profile_id,
#       sim_df = list(sim_df),
#       time_sd = time_sd,
#       mean_mel_sd = mean_mel_sd,
#       min_mel_sd = min_mel_sd,
#       max_mel_sd = max_mel_sd,
#       dlmo_est = result$ip$inflection_point_fine$x,
#       dlmo_full_result = list(result),
#       condition = condition,
#       replicate = replicate
#     )
#   }
#
#   conditions <- list(
#     clean = list(time_sd = 0, add_mel = FALSE),
#     mel_only = list(time_sd = 0, add_mel = TRUE),
#     time_5min = list(time_sd = time_sd_5, add_mel = FALSE),
#     time_10min = list(time_sd = time_sd_10, add_mel = FALSE),
#     time_20min = list(time_sd = time_sd_20, add_mel = FALSE),
#     both_axes_10 = list(time_sd = time_sd_10, add_mel = TRUE)
#   )
#
#   all_results <- purrr::map_dfr(names(conditions), function(cond) {
#     p <- conditions[[cond]]
#     reps <- if (cond == "clean") 1 else n_rep
#     purrr::map_dfr(seq_len(reps), function(i) {
#       run_one(p$time_sd, p$add_mel, cond, i)
#     })
#   })
#
#   list(results = all_results, errors = dplyr::bind_rows(error_rows))
# }
#
# library(dplyr)
# library(purrr)
# library(dlmoR)
# library(readr)
#
# # ────────────────────────────────────────────────────
# # SET PROFILE PATH AND ID
# # ────────────────────────────────────────────────────
# profile_path <- "inst/extdata/first_Visit_1_PB_8.csv"   # <-- Replace with your actual path
# profile_id   <- "first_Visit_1_PB_8"                    # <-- Base name, no extension
#
# # ────────────────────────────────────────────────for (i in seq(2, length(times_dec))) {
# # ────────────────────────────────────────────────────
# ────
# # RUN THE PROCESSING FUNCTION
# # ────────────────────────────────────────────────────
# result_fix <- process_single_profile(profile_id, profile_path, n_rep = 10)
#
# # ────────────────────────────────────────────────────
# # SAVE THE RESULT
# # ────────────────────────────────────────────────────
# out_dir <- "outputs/noisy_dlmo_results"
# dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
#
# saveRDS(result_fix, file.path(out_dir, paste0(profile_id, ".rds")))
#
# iteration_errors_df <- result_fix$errors
# readr::write_csv(iteration_errors_df, file.path(out_dir, "noise_simulation_iteration_errors.csv"))
#
# # Save errors (including sim_df!)
# readr::write_rds(result_fix$errors, file.path(out_dir, paste0(profile_id, "_iteration_errors.rds")))
#
# # OR save as CSV if you drop sim_df (list-column not supported in CSV)
# readr::write_csv(dplyr::select(result_fix$errors, -sim_df), file.path(out_dir, paste0(profile_id, "_iteration_errors.csv")))
#
#
# # (Optional) view a summary
# # print(result)
#
