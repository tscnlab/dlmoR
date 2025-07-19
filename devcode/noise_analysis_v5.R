library(dplyr)
library(purrr)
library(dlmoR)
library(lubridate)
library(hms)

# --- Determine which profiles to process -----------------------------------
profile_folder <- system.file("extdata/", package = "dlmoR")
csv_files      <- list.files(profile_folder, pattern = "\\.csv$", full.names = TRUE)

# Results directory (ensure it matches the one you save into)
out_dir        <- "outputs/noisy_dlmo_results"
result_files   <- list.files(out_dir, pattern = "\\.rds$")
processed_ids  <- tools::file_path_sans_ext(basename(result_files))

# Named vector of all profiles, filter out already processed
profiles <- set_names(csv_files,
                      tools::file_path_sans_ext(basename(csv_files)))
profiles <- profiles[!names(profiles) %in% processed_ids]

# --- Helper: enforce minimum gap on decimal times --------------------------
enforce_min_gap <- function(times_dec, min_gap = 1/60) {
  for (i in seq(2, length(times_dec))) {
    gap <- times_dec[i] - times_dec[i - 1]
    if (gap < min_gap) times_dec[i] <- times_dec[i - 1] + min_gap
  }
  times_dec
}

# --- Run one profile with full outputs --------------------------------------
run_profile_noise_full <- function(sample_dlmo,
                                   n_rep          = 20,
                                   mel_sd         = 0.3,
                                   time_sd_main   = NULL,
                                   extra_time_sds = NULL,
                                   seed           = 42) {
  df        <- sample_dlmo$prof
  ref_times <- df$datetime
  dec_t     <- posixct_to_decimal(ref_times, ref_times)

  # Determine time SDs if not provided
  min_int <- min(diff(sort(dec_t)))
  if (is.null(time_sd_main)) time_sd_main <- min_int / 2
  if (is.null(extra_time_sds)) extra_time_sds <- c(time_sd_main/2, time_sd_main*2)

  # Interpolation
  interp_fun <- splinefun(x = dec_t, y = df$melatonin, method = "monoH.FC")
  set.seed(seed)

  # Runner for one perturbed replicate
  do_one_full <- function(time_sd = 0, mel_sd = 0) {
    t_pert <- dec_t + rnorm(length(dec_t), 0, time_sd)
    t_pert <- enforce_min_gap(sort(t_pert), min_gap = 1/60)

    mel_int  <- interp_fun(t_pert)
    mel_pert <- mel_int + rnorm(length(mel_int), 0, mel_sd)

    sim_df <- tibble(
      datetime  = decimal_to_posixct(t_pert, ref_times),
      melatonin = mel_pert
    )

    tryCatch(calculate_dlmo(sim_df), error = function(e) NULL)
  }

  # Clean output
  clean_out <- calculate_dlmo(df)

  # Monte Carlo lists
  mel_list  <- replicate(n_rep, do_one_full(time_sd = 0, mel_sd = mel_sd), simplify = FALSE)
  time_sds  <- c(time_sd_main, extra_time_sds)
  time_list <- map(time_sds, ~ replicate(n_rep,
                                         do_one_full(time_sd = .x, mel_sd = 0), simplify = FALSE))
  both_list <- replicate(n_rep, do_one_full(time_sd = time_sd_main, mel_sd = mel_sd), simplify = FALSE)

  # Assemble results
  types <- c(
    "clean",
    rep("mel_only",      n_rep),
    rep(paste0("time_sd=", time_sds[1]), n_rep),
    rep(paste0("time_sd=", time_sds[2]), n_rep),
    rep(paste0("time_sd=", time_sds[3]), n_rep),
    rep("both_axes",      n_rep)
  )

  time_sds_col <- c(0,
                    rep(0,                n_rep),
                    rep(time_sds[1],      n_rep),
                    rep(time_sds[2],      n_rep),
                    rep(time_sds[3],      n_rep),
                    rep(time_sd_main,     n_rep)
  )

  mel_sds_col <- c(0,
                   rep(mel_sd,           n_rep),
                   rep(0,                n_rep*3),
                   rep(mel_sd,           n_rep)
  )

  full_outputs <- c(
    list(clean_out),
    mel_list,
    time_list[[1]],
    time_list[[2]],
    time_list[[3]],
    both_list
  )

  tibble(
    profile_id  = sample_dlmo$profile_name,
    type        = types,
    time_sd     = time_sds_col,
    mel_sd      = mel_sds_col,
    full_output = full_outputs
  ) %>%
    mutate(
      dlmo_est = map_dbl(full_output, ~ .x$ip$inflection_point_fine$x %||% NA_real_)
    )
}

# --- Build sample_dlmo list from CSVs ---------------------------------------
all_samples <- imap(profiles, function(path, profile_name) {
  raw_df <- read.csv(path) %>%
    mutate(
      datetime  = as.POSIXct(datetime, tz = "UTC"),
      melatonin = melatonin
    )

  sample <- calculate_dlmo(raw_df)
  sample$prof         <- raw_df
  sample$profile_name <- profile_name
  sample
})

# --- Run over all and save ----------------------------------------------

# To process only a subset of profiles, first select their names or indices. For example, to run on profiles "A" and "B":
# subset_samples <- all_samples[c("A", "B")]
# Or by index: subset_samples <- all_samples[1:3]
# Then call with a custom number of Monte Carlo repetitions:
# run_all_profiles_full(subset_samples, out_dir = out_dir, n_rep = 50)

# To run on all profiles with default settings:
run_all_profiles_full(all_samples[1:2], out_dir = out_dir, n_rep = 2)
