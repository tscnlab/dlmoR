library(dplyr)
library(purrr)
library(readr)

# Set input and output directories
input_dir  <- "~/Documents/Projects/DLMO/dlmoRpaperresults/Blume/noise_results/"

output_dir <- "~/Documents/Projects/DLMO/dlmoRpaperresults/Blume/noise_results_skinny/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# List all large .rds files
rds_files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)

# Loop through files and clean
walk(rds_files, function(file) {
  message("Cleaning: ", basename(file))

  # Try to read and strip the large column
  safe_data <- tryCatch(readRDS(file), error = function(e) {
    message("Failed to read: ", file, " (", e$message, ")")
    return(NULL)
  })

  if (!is.null(safe_data) && "dlmo_full_result" %in% names(safe_data)) {
    safe_data <- select(safe_data, -dlmo_full_result)
  }

  # Save cleaned version to new directory
  new_name <- file.path(output_dir, basename(file))
  saveRDS(safe_data, new_name)
})
