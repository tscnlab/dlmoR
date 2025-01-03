library(dplyr)
library(readr)
library(tidyr)

# Read the data from CSV
filename <- 'C:/Users/FT/Documents/Salma/tSCN/Projects/dlmoR_nonpackage_files/BlumeetalData/Data/melatonin_data.csv'  # Replace with the path to your CSV file
data <- read.csv(filename)

# Ensure PB_No and Visit columns are properly recognized
data$PB_No <- as.factor(data$PB_No)
data$Visit <- as.factor(data$Visit)

# Ensure datetime sequence
start_time <- as.POSIXct("2025-01-01 00:00:00")
datetime_seq <- seq(from = start_time, by = "30 min", length.out = 14)

# Create a function to split and save the data
split_and_save <- function(df) {
  # Check if PB_No and Visit columns exist
  if (!("PB_No" %in% colnames(df)) || !("Visit" %in% colnames(df))) {
    stop("PB_No or Visit columns are missing from the data.")
  }

  # Split by Visit and PB_No
  split_data <- split(df, list(df$PB_No, df$Visit), drop = TRUE)

  # Iterate over each subset
  for (key in names(split_data)) {
    subset <- split_data[[key]]

    # Sort by Number (just in case)
    subset <- subset[order(subset$Number), ]

    # Generate filenames
    base_filename <- paste0("Visit_", subset$Visit[1], "_PB_", subset$PB_No[1], ".csv")
    dir.create("output", showWarnings = FALSE)

    # Split into two parts (first 14 rows and remaining)
    if (nrow(subset) >= 28) {
      first_part <- subset[1:14, ]
      second_part <- subset[(nrow(subset) - 13):nrow(subset), ]

      # Apply datetime explicitly
      first_part$datetime <- datetime_seq
      second_part$datetime <- datetime_seq

      # Format as dataframe
      first_part <- data.frame(datetime = format(first_part$datetime, "%Y-%m-%d %H:%M:%S"), melatonin = first_part$Value)
      second_part <- data.frame(datetime = format(second_part$datetime, "%Y-%m-%d %H:%M:%S"), melatonin = second_part$Value)

      # Save the files with semicolon delimiter and suppress quotes
      write.table(first_part, file = paste0("output/first_", base_filename), row.names = FALSE, sep = ";", dec = ".", quote = FALSE, col.names = TRUE)
      write.table(second_part, file = paste0("output/second_", base_filename), row.names = FALSE, sep = ";", dec = ".", quote = FALSE, col.names = TRUE)
    } else {
      warning(paste("Skipping", base_filename, "- not enough rows."))
    }
  }
}

# Run the split function
split_and_save(data)

# # Process the output files for DLMO calculation
# calculate_dlmo_for_files <- function(output_dir = "output", threshold = 5) {
#   files <- list.files(output_dir, pattern = "*.csv", full.names = TRUE)
#
#   for (file in files) {
#     file_name <- tools::file_path_sans_ext(basename(file))
#     dlmo_result <- calculate_dlmo(file_path = file, threshold = threshold)
#     assign(paste0(file_name, "_dlmo"), dlmo_result, envir = .GlobalEnv)
#   }
# }
#
# # Run DLMO calculation
# calculate_dlmo_for_files()


# # Process the output files for DLMO calculation and save plots
# calculate_dlmo_for_files <- function(output_dir = "output", threshold = 5) {
#   files <- list.files(output_dir, pattern = "*.csv", full.names = TRUE)
#
#   for (file in files) {
#     file_name <- tools::file_path_sans_ext(basename(file))
#     dlmo_result <- calculate_dlmo(file_path = file, threshold = threshold)
#     assign(paste0(file_name, "_dlmo"), dlmo_result, envir = .GlobalEnv)
#
#     # Save coarse plot as PDF
#     pdf(paste0("output/plots", file_name, "_dlmo_coarse_plot.pdf"))
#     print(dlmo_result$dlmoplotcoarse)
#     dev.off()
#
#     # Save fine plot as PDF
#     pdf(paste0("output/plots", file_name, "_dlmo_fine_plot.pdf"))
#     print(dlmo_result$dlmoplotfine)
#     dev.off()
#   }
# }
#
# # Run DLMO calculation
# calculate_dlmo_for_files()




# calculate_dlmo_for_files <- function(output_dir = "output", threshold = 5) {
#   files <- list.files(output_dir, pattern = "*.csv", full.names = TRUE)
#
#   # Process only the first 2 files
#   selected_files <- files[1:2]
#
#   for (file in selected_files) {
#     file_name <- tools::file_path_sans_ext(basename(file))
#     dlmo_result <- calculate_dlmo(file_path = file, threshold = threshold)
#     assign(paste0(file_name, "_dlmo"), dlmo_result, envir = .GlobalEnv)
#
#     # Save coarse plot as PDF
#     pdf(paste0("output/plots/", file_name, "_dlmo_coarse_plot.pdf"))
#     print(dlmo_result$dlmoplotcoarse)
#     dev.off()
#
#     # Save fine plot as PDF
#     pdf(paste0("output/plots/", file_name, "_dlmo_fine_plot.pdf"))
#     print(dlmo_result$dlmoplotfine)
#     dev.off()
#   }
# }
#
# # Run DLMO calculation
# calculate_dlmo_for_files()


# calculate_dlmo_for_files <- function(output_dir = "output", threshold = 5) {
#   files <- list.files(output_dir, pattern = "*.csv", full.names = TRUE)
#
#   # Create the plots subdirectory if it doesn't exist
#   dir.create(file.path(output_dir, "plots"), showWarnings = FALSE)
#
#   # Process only the first 2 files
#   #selected_files <- files[1:2]
#
#   for (file in files) {
#     file_name <- tools::file_path_sans_ext(basename(file))
#     dlmo_result <- calculate_dlmo(file_path = file, threshold = threshold)
#     assign(paste0(file_name, "_dlmo"), dlmo_result, envir = .GlobalEnv)
#
#     # Save coarse plot as PDF
#     pdf(paste0(output_dir, "/plots/", file_name, "_dlmo_coarse_plot.pdf"))
#     print(dlmo_result$dlmoplotcoarse)
#     dev.off()
#
#     # Save fine plot as PDF
#     pdf(paste0(output_dir, "/plots/", file_name, "_dlmo_fine_plot.pdf"))
#     print(dlmo_result$dlmoplotfine)
#     dev.off()
#
#     # Store DLMO inflection point
#     results <- rbind(results, data.frame(data_source = file_name, dlmoR_ip_decimalhours = dlmo_result$ip$inflection_point$x))
#   }
#   # Save results to CSV
#   write.csv(results, file = "output/dlmo_results_summary.csv", row.names = FALSE)
# }
#
#
# calculate_dlmo_for_files()

########## THIS WORKS ##################
# Initialize results dataframe globally
# results <- data.frame(data_source = character(), dlmoR_ip_decimalhours = numeric(), stringsAsFactors = FALSE)
#
# calculate_dlmo_for_files <- function(output_dir = "output", threshold = 5) {
#   files <- list.files(output_dir, pattern = "*.csv", full.names = TRUE)
# }
#   # Create the plots subdirectory if it doesn't exist
#   dir.create(file.path(output_dir, "plots"), showWarnings = FALSE)
#
#   for (file in files) {
#     file_name <- tools::file_path_sans_ext(basename(file))
#     dlmo_result <- calculate_dlmo(file_path = file, threshold = threshold)
#     assign(paste0(file_name, "_dlmo"), dlmo_result, envir = .GlobalEnv)
#
#     # Save coarse plot as PDF
#     pdf(paste0(output_dir, "/plots/", file_name, "_dlmo_coarse_plot.pdf"))
#     print(dlmo_result$dlmoplotcoarse)
#     dev.off()
#
#     # Save fine plot as PDF
#     pdf(paste0(output_dir, "/plots/", file_name, "_dlmo_fine_plot.pdf"))
#     print(dlmo_result$dlmoplotfine)
#     dev.off()
#
#     # Store DLMO inflection point
#     new_row <- data.frame(data_source = file_name, dlmoR_ip_decimalhours = dlmo_result$ip$inflection_point$x)
#     assign("results", rbind(get("results", envir = .GlobalEnv), new_row), envir = .GlobalEnv)
#   }
#
#   # Save results to CSV
#   write.csv(results, file = "output/dlmo_results_summary.csv", row.names = FALSE)
# }
#
# # Run DLMO calculation
# calculate_dlmo_for_files()


########### this allows error handling
results <- data.frame(data_source = character(), dlmoR_ip_decimalhours = numeric(), stringsAsFactors = FALSE)
error_files <- list()  # List to store files that encounter errors

calculate_dlmo_for_files <- function(output_dir = "output", threshold = 5) {
  files <- list.files(output_dir, pattern = "*.csv", full.names = TRUE)

  # Create the plots subdirectory if it doesn't exist
  dir.create(file.path(output_dir, "plots"), showWarnings = FALSE)

  for (file in files) {
    file_name <- tools::file_path_sans_ext(basename(file))

    # Wrap the processing inside tryCatch to handle errors
    tryCatch({
      dlmo_result <- calculate_dlmo(file_path = file, threshold = threshold)
      assign(paste0(file_name, "_dlmo"), dlmo_result, envir = .GlobalEnv)

      # Save coarse plot as PDF
      pdf(paste0(output_dir, "/plots/", file_name, "_dlmo_coarse_plot.pdf"))
      print(dlmo_result$dlmoplotcoarse)
      dev.off()

      # Save fine plot as PDF
      pdf(paste0(output_dir, "/plots/", file_name, "_dlmo_fine_plot.pdf"))
      print(dlmo_result$dlmoplotfine)
      dev.off()

      # Store DLMO inflection point
      new_row <- data.frame(data_source = file_name, dlmoR_ip_decimalhours = dlmo_result$ip$inflection_point$x)
      assign("results", rbind(get("results", envir = .GlobalEnv), new_row), envir = .GlobalEnv)

    }, error = function(e) {
      # In case of an error, print the file name and the error message
      print(paste("Error occurred with file:", file_name))
      print(e$message)

      # Optionally, store the file name in a list for later review
      error_files <<- c(error_files, file_name)
    })
  }

  # Save results to CSV
  write.csv(results, file = "output/dlmo_results_summary.csv", row.names = FALSE)

  # If there were any errors, print them out
  if (length(error_files) > 0) {
    cat("Errors occurred with the following files:\n")
    print(error_files)
  } else {
    cat("All files processed successfully.\n")
  }
}

# Run DLMO calculation
calculate_dlmo_for_files()

