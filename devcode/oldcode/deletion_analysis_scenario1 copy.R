# -----------------------------
# Required packages
# -----------------------------
library(readr)
library(dplyr)
library(purrr)
library(ggplot2)

# -----------------------------
# Helper functions
# -----------------------------
extract_dlmo_value <- function(dlmo_result) {
  dlmo_result$ip$inflection_point_fine$x
}

relative_minutes_to_dlmo <- function(timestamps, dlmo_time) {
  as.numeric(difftime(timestamps, dlmo_time, units = "mins"))
}

# -----------------------------
# Load a single profile
# -----------------------------
profile_folder <- system.file("extdata/", package = "dlmoR")  # <- replace if needed
csv_files <- list.files(profile_folder, pattern = "\\.csv$", full.names = TRUE)

df <- read_csv(csv_files[1], show_col_types = FALSE)
profile_id <- tools::file_path_sans_ext(basename(csv_files[1]))

# -----------------------------
# Run Scenario 1 on this profile
# -----------------------------
full_dlmo_result <- tryCatch({
  calculate_dlmo(df, threshold = 5)
}, error = function(e) {
  message("Full DLMO failed: ", e$message)
  return(NULL)
})

if (is.null(full_dlmo_result)) {
  message("Skipping profile ", profile_id, " due to DLMO error on full profile.")
} else {
  full_dlmo <- extract_dlmo_value(full_dlmo_result)
  full_dlmo_time <- full_dlmo_result$ip$inflection_point_fine$datetime

  scenario1 <- map(seq_len(nrow(df)), function(i) {
    df_deleted <- df[-i, ]
    cat("Testing deletion at index:", i, "\n")

    tryCatch({
      dlmo_deleted <- extract_dlmo_value(calculate_dlmo(df_deleted, threshold = 5))
      delta <- dlmo_deleted - full_dlmo
      cat(sprintf("Success at index %d: Δ DLMO = %.4f\n", i, delta))
      # print(str(list(
      #   profile = profile_id,
      #   deleted_time = df$datetime[i],
      #   delta_minutes_from_dlmo = list(relative_minutes_to_dlmo(df$datetime[i], full_dlmo_time)),
      #   delta_dlmo = delta,
      #   dlmo_success = TRUE
      # )))
      out <- tibble(
        profile = profile_id,
        deleted_time = df$datetime[i],
        delta_minutes_from_dlmo = relative_minutes_to_dlmo(df$datetime[i], full_dlmo_time),
        delta_dlmo = delta,
        dlmo_success = TRUE
      )
      print(out)
      return(out)

    }, error = function(e) {
      cat(sprintf("Failure at index %d: %s\n", i, e$message))
      tibble(
        profile = profile_id,
        deleted_time = df$datetime[i],
        delta_minutes_from_dlmo = relative_minutes_to_dlmo(df$datetime[i], full_dlmo_time),
        delta_dlmo = NA_real_,
        dlmo_success = FALSE
      )
    })
  }) %>% bind_rows()



  print(scenario1)

  # Optional: Plot
  ggplot(scenario1, aes(x = delta_minutes_from_dlmo, y = delta_dlmo)) +
    geom_point(aes(color = dlmo_success)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
      title = paste("Δ DLMO after single-point deletions:", profile_id),
      x = "Time of deleted point (minutes relative to full-profile DLMO)",
      y = "Change in DLMO (hours)"
    ) +
    theme_minimal()
}



# Minimal diagnostic version of Scenario 1 for a single profile
run_single_point_diagnostics <- function(profile_id, df) {
  full_dlmo_result <- tryCatch({
    calculate_dlmo(df, threshold = 5)
  }, error = function(e) {
    message("Full DLMO failed: ", e$message)
    return(NULL)
  })

  if (is.null(full_dlmo_result)) {
    message("Skipping profile ", profile_id, " due to DLMO error on full profile.")
    return(tibble())
  }

  full_dlmo <- extract_dlmo_value(full_dlmo_result)
  full_dlmo_time <- full_dlmo_result$ip$inflection_point_fine$datetime

  scenario1_list <- list()

  for (i in seq_len(nrow(df))) {
    cat("Testing deletion at index:", i, "\n")
    df_deleted <- df[-i, ]

    result <- tryCatch({
      dlmo_deleted <- extract_dlmo_value(calculate_dlmo(df_deleted, threshold = 5))
      delta <- dlmo_deleted - full_dlmo
      row <- tibble(
        profile = profile_id,
        deleted_time = df$datetime[i],
        delta_minutes_from_dlmo = relative_minutes_to_dlmo(df$datetime[i], full_dlmo_time),
        delta_dlmo = delta,
        dlmo_success = TRUE
      )
      print(row)
      row
    }, error = function(e) {
      message("Failure at index ", i, ": ", e$message)
      tibble(
        profile = profile_id,
        deleted_time = df$datetime[i],
        delta_minutes_from_dlmo = relative_minutes_to_dlmo(df$datetime[i], full_dlmo_time),
        delta_dlmo = NA_real_,
        dlmo_success = FALSE
      )
    })

    scenario1_list[[i]] <- result
  }

  scenario1 <- bind_rows(scenario1_list)
  return(scenario1)
}

