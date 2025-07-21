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

  clean_res <- calculate_dlmo(df, threshold = 5)
  clean_time <- clean_res$ip$inflection_point_fine$x
  dec_t <- posixct_to_decimal(df$datetime, df$datetime)
  interp_fun <- splinefun(dec_t, df$melatonin, method = "monoH.FC")

  time_sd_5  <- 5 / 60
  time_sd_10 <- 10 / 60
  time_sd_20 <- 20 / 60

  run_one <- function(time_sd, add_mel_noise) {
    t_pert <- enforce_min_gap(sort(dec_t + rnorm(length(dec_t), 0, time_sd)))
    mel_int <- interp_fun(t_pert)

    if (add_mel_noise) {
      sd_vec <- (cv / 100) * mel_int
      mel_pert <- mel_int + rnorm(length(mel_int), 0, sd_vec)
      mean_mel_sd <- mean(sd_vec)
      min_mel_sd <- min(sd_vec)
      max_mel_sd <- max(sd_vec)
    } else {
      mel_pert <- mel_int
      mean_mel_sd <- min_mel_sd <- max_mel_sd <- 0
    }

    sim_df <- tibble(
      datetime  = decimal_to_posixct(t_pert, df$datetime),
      melatonin = mel_pert
    )

    result <- tryCatch(calculate_dlmo(sim_df, threshold = 5), error = function(e) NULL)

    tibble(
      profile = profile_id,
      time_sd = time_sd,
      mean_mel_sd = mean_mel_sd,
      min_mel_sd = min_mel_sd,
      max_mel_sd = max_mel_sd,
      dlmo_est = ifelse(is.null(result), NA_real_, result$ip$inflection_point_fine$x),
      dlmo_full_result = list(result)
    )
  }

  conditions <- list(
    clean = list(time_sd = 0, add_mel = FALSE),
    mel_only = list(time_sd = 0, add_mel = TRUE),
    time_5min = list(time_sd = time_sd_5, add_mel = FALSE),
    time_10min = list(time_sd = time_sd_10, add_mel = FALSE),
    time_20min = list(time_sd = time_sd_20, add_mel = FALSE),
    both_axes_10 = list(time_sd = time_sd_10, add_mel = TRUE)
  )

  map_dfr(names(conditions), function(cond) {
    p <- conditions[[cond]]
    reps <- if (cond == "clean") 1 else n_rep
    map_dfr(seq_len(reps), function(i) {
      run_one(p$time_sd, p$add_mel) %>%
        mutate(condition = cond, replicate = i)
    })
  })
}

# ────────────────────────────────────────────────────
# RUN IN PARALLEL (with informative printouts)
# ────────────────────────────────────────────────────
safe_process <- safely(process_single_profile)

with_progress({
  p <- progressor(along = profiles_to_run)

  results <- future_imap(profiles_to_run, function(path, profile_id) {
    res <- safe_process(profile_id, path)

    if (!is.null(res$result)) {
      saveRDS(res$result, file.path(out_dir, paste0(profile_id, ".rds")))
      message("[SUCCESS] Completed profile: ", profile_id)
    } else {
      message("[ERROR] Failed profile: ", profile_id, " - ", res$error$message)
    }

    p()
    list(profile = profile_id, error = res$error)
  }, .options = furrr_options(seed = TRUE))
})

# ────────────────────────────────────────────────────
# ERROR LOG CSV
# ────────────────────────────────────────────────────
error_log <- keep(results, ~ !is.null(.x$error)) %>%
  map_dfr(~ tibble(profile = .x$profile, error_message = .x$error$message))

write_csv(error_log, file.path(out_dir, "noise_simulation_errors.csv"))

# ────────────────────────────────────────────────────
# COMBINED RESULTS
# ────────────────────────────────────────────────────
# all_result_files <- list.files(out_dir, pattern = "\\.rds$", full.names = TRUE)
# all_results <- map_dfr(all_result_files, readRDS)

# saveRDS(all_results, file.path(out_dir, "dlmo_noise_all_results.rds"))
