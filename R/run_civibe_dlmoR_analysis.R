# script for running civibe data
# results <- data.frame(data_source = character(), dlmoR_ip_decimalhours = numeric(), stringsAsFactors = FALSE)
# error_files <- list()  # List to store files that encounter errors
#
# calculate_dlmo_for_files <- function(output_dir = "output/civibe", threshold = 5) {
#   files <- list.files(output_dir, pattern = "*.csv", full.names = TRUE)
#
#   # Create the plots subdirectory if it doesn't exist
#   dir.create(file.path(output_dir, "plots"), showWarnings = FALSE)
#
#   for (file in files) {
#     file_name <- tools::file_path_sans_ext(basename(file))
#
#     # Wrap the processing inside tryCatch to handle errors
#     tryCatch({
#       dlmo_result <- calculate_dlmo(file_path = file, threshold = threshold)
#       assign(paste0(file_name, "_dlmo"), dlmo_result, envir = .GlobalEnv)
#
#       # Save coarse plot as PDF
#       pdf(paste0(output_dir, "/plots/", file_name, "_dlmo_coarse_plot.pdf"))
#       print(dlmo_result$dlmoplotcoarse)
#       dev.off()
#
#       # Save fine plot as PDF
#       pdf(paste0(output_dir, "/plots/", file_name, "_dlmo_fine_plot.pdf"))
#       print(dlmo_result$dlmoplotfine)
#       dev.off()
#
#       # Store DLMO inflection point
#       new_row <- data.frame(data_source = file_name, dlmoR_ip_decimalhours = dlmo_result$ip$inflection_point$x)
#       assign("results", rbind(get("results", envir = .GlobalEnv), new_row), envir = .GlobalEnv)
#
#     }, error = function(e) {
#       # In case of an error, print the file name and the error message
#       print(paste("Error occurred with file:", file_name))
#       print(e$message)
#
#       # Optionally, store the file name in a list for later review
#       error_files <<- c(error_files, file_name)
#     })
#   }
#
#   # Save results to CSV
#   write.csv(results, file = "output/civibe/civibe_dlmo_results_summary.csv", row.names = FALSE)
#
#   # If there were any errors, print them out
#   if (length(error_files) > 0) {
#     cat("Errors occurred with the following files:\n")
#     print(error_files)
#   } else {
#     cat("All files processed successfully.\n")
#   }
# }
#
# # Run DLMO calculation
# calculate_dlmo_for_files()


# Script for running civibe data from a list
# results <- data.frame(data_source = character(), dlmoR_ip_decimalhours = numeric(), stringsAsFactors = FALSE)
# error_files <- list()  # List to store files that encounter errors
#
# calculate_dlmo_for_list <- function(data_list, output_dir = "output/civibe", threshold = 2.3) {
#   # Create the output and plots directories if they don't exist
#   dir.create(output_dir, showWarnings = FALSE)
#   dir.create(file.path(output_dir, "plots"), showWarnings = FALSE)
#
#   for (data_name in names(data_list)) {
#     tryCatch({
#       # Run calculate_dlmo directly on the list element
#       dlmo_result <- calculate_dlmo(data = data_list[[data_name]], threshold = threshold)
#       assign(paste0(data_name, "_dlmo"), dlmo_result, envir = .GlobalEnv)
#       print(dlmo_result)
#       # Save coarse plot as PDF
#       pdf(paste0(output_dir, "/plots/", data_name, "_dlmo_coarse_plot.pdf"))
#       print(dlmo_result$dlmoplotcoarse)
#       dev.off()
#
#       # Save fine plot as PDF
#       pdf(paste0(output_dir, "/plots/", data_name, "_dlmo_fine_plot.pdf"))
#       print(dlmo_result$dlmoplotfine)
#       dev.off()
#
#       # Store DLMO inflection point
#       new_row <- data.frame(data_source = data_name, dlmoR_ip_decimalhours = dlmo_result$ip$inflection_point$x)
#       assign("results", rbind(get("results", envir = .GlobalEnv), new_row), envir = .GlobalEnv)
#
#     }, error = function(e) {
#       # In case of an error, print the data name and the error message
#       print(paste("Error occurred with data:", data_name))
#       print(e$message)
#
#       # Store the data name in a list for later review
#       error_files <<- c(error_files, data_name)
#     })
#   }
#
#   # Save results to CSV
#   write.csv(results, file = paste0(output_dir, "/civibe_dlmo_results_summary.csv"), row.names = FALSE)
#
#   # If there were any errors, print them out
#   if (length(error_files) > 0) {
#     cat("Errors occurred with the following data:\n")
#     print(error_files)
#   } else {
#     cat("All data processed successfully.\n")
#   }
# }
#
# # Run DLMO calculation
# calculate_dlmo_for_list(civibe_melatonin_FD_dh_daysplit, threshold = 2.3)
#

# results <- data.frame(data_source = character(), dlmoR_ip_decimalhours = numeric(), stringsAsFactors = FALSE)
# error_files <- list()  # List to store data that encounter errors
#
# calculate_dlmo_for_list <- function(data_list, output_dir = "output/civibe", threshold = 2.3) {
#   # Create the output and plots directories if they don't exist
#   dir.create(output_dir, showWarnings = FALSE)
#   dir.create(file.path(output_dir, "plots"), showWarnings = FALSE)
#
#   for (data_name in names(data_list)) {
#     temp_data <- data_list[[data_name]]  # Assign list element to a temp variable
#
#     # Wrap the processing inside tryCatch to handle errors
#     tryCatch({
#       # Run calculate_dlmo directly on the temp data
#       dlmo_result <- calculate_dlmo(data = temp_data, threshold = threshold)
#       assign(paste0(data_name, "_dlmo"), dlmo_result, envir = .GlobalEnv)
#
#       # Save coarse plot as PDF
#       pdf(paste0(output_dir, "/plots/", data_name, "_dlmo_coarse_plot.pdf"))
#       print(dlmo_result$dlmoplotcoarse)
#       dev.off()
#
#       # Save fine plot as PDF
#       pdf(paste0(output_dir, "/plots/", data_name, "_dlmo_fine_plot.pdf"))
#       print(dlmo_result$dlmoplotfine)
#       dev.off()
#
#       # Store DLMO inflection point
#       new_row <- data.frame(data_source = data_name, dlmoR_ip_decimalhours = dlmo_result$ip$inflection_point$x)
#       assign("results", rbind(get("results", envir = .GlobalEnv), new_row), envir = .GlobalEnv)
#
#     }, error = function(e) {
#       # In case of an error, print the data name and the error message
#       print(paste("Error occurred with data:", data_name))
#       print(e$message)
#
#       # Store the data name in a list for later review
#       error_files <<- c(error_files, data_name)
#     })
#   }
#
#   # Save results to CSV
#   write.csv(results, file = paste0(output_dir, "/civibe_dlmo_results_summary.csv"), row.names = FALSE)
#
#   # If there were any errors, print them out
#   if (length(error_files) > 0) {
#     cat("Errors occurred with the following data:\n")
#     print(error_files)
#   } else {
#     cat("All data processed successfully.\n")
#   }
# }
#
# # Run DLMO calculation
# calculate_dlmo_for_list(civibe_melatonin_FD_dh_daysplit, threshold = 2.3)
#

results <- data.frame(
  data_source = character(),
  dlmoR_ip_decimalhours = numeric(),
  stringsAsFactors = FALSE
)
error_files <- list()  # List to store data that encounter errors

calculate_dlmo_for_list <- function(data_list, output_dir = "output/civibe", threshold = 2.3) {
  # Create the output and plots directories if they don't exist
  dir.create(output_dir, showWarnings = FALSE)
  dir.create(file.path(output_dir, "plots"), showWarnings = FALSE)

  for (data_name in names(data_list)) {
    temp_data <- data_list[[data_name]]  # Assign list element to a temp variable
    tryCatch({
      # Run calculate_dlmo directly on the temp data
      dlmo_result <- calculate_dlmo(data = temp_data, threshold = threshold)
      # print("dlmoresult")
      # print(dlmo_result)
      assign(paste0(data_name, "_dlmo"), dlmo_result, envir = .GlobalEnv)

      # Save coarse plot as PDF
      pdf(paste0(output_dir, "/plots/", data_name, "_dlmo_coarse_plot.pdf"))
      print(dlmo_result$dlmoplotcoarse)
      dev.off()

      # Save fine plot as PDF
      pdf(paste0(output_dir, "/plots/", data_name, "_dlmo_fine_plot.pdf"))
      print(dlmo_result$dlmoplotfine)
      dev.off()

      # Extract the inflection point
      inflection_value <- dlmo_result$ip$inflection_point$x

      # Debugging print
      print(paste("Processing:", data_name))
      print(paste("Inflection point (x):", inflection_value))

      # Check if the inflection point is valid
      if (!is.null(inflection_value) && is.numeric(inflection_value)) {
        new_row <- data.frame(
          data_source = as.character(data_name),
          dlmoR_ip_decimalhours = as.numeric(inflection_value),
          stringsAsFactors = FALSE
        )

        # Append with tryCatch to avoid halting on errors
        tryCatch({
          results <<- rbind(results, new_row)
        }, error = function(e) {
          print(paste("Failed to append data for:", data_name))
          print(e$message)
          error_files <<- c(error_files, data_name)
        })
      } else {
        stop(paste("Inflection point missing or invalid for:", data_name))
      }

    }, error = function(e) {
      # Print the data name and error message
      print(paste("Error occurred with data:", data_name))
      print(e$message)

      # Store the data name in a list for later review
      error_files <<- c(error_files, data_name)
    })
  }

  # Save results to CSV
  write.csv(results, file = paste0(output_dir, "/civibe_dlmo_results_summary.csv"), row.names = FALSE)

  # If there were any errors, print them out
  if (length(error_files) > 0) {
    cat("Errors occurred with the following data:\n")
    print(error_files)
  } else {
    cat("All data processed successfully.\n")
  }
}

# Run DLMO calculation
calculate_dlmo_for_list(civibe_melatonin_FD_dh_daysplit, threshold = 2.3)
