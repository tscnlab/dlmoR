define_intermediate_segment <- function(profile_data, threshold = threshold) {
  # Ensure the data is sorted by datetime
  profile_data <- profile_data %>% dplyr::arrange(.data$datetime)

  # Find the last row of the base part
  last_base_row <- max(which(profile_data$base == 1), na.rm = TRUE)

  # Find the first row of the ascending part
  first_ascending_row <- min(which(profile_data$ascending == 1), na.rm = TRUE)

  # Check if there are any intermediate rows
  if (last_base_row < first_ascending_row - 1) {
    # Identify intermediate rows (between last base and first ascending)
    intermediate_rows <- (last_base_row + 1):(first_ascending_row - 1)

    # Create a new column 'intermediate', defaulting to 0
    profile_data <- profile_data %>%
      dplyr::mutate(intermediate = 0)

    # Assign '1' to the intermediate rows
    profile_data$intermediate[intermediate_rows] <- 1
  }

  # Return the updated profile data
  return(profile_data)
}


#### THIS BELOW FAILS EXCEPT IN CASE OF NO BASE!!!!! FIX IT #############

# define_intermediate_segment <- function(profile_data, threshold) {
#   # Ensure the data is sorted by datetime
#   profile_data <- profile_data %>% dplyr::arrange(.data$datetime)
#
#   # Check if there are any base or ascending points
#   base_exists <- any(profile_data$base == 1, na.rm = TRUE)
#   ascending_exists <- any(profile_data$ascending == 1, na.rm = TRUE)
#
#   if (ascending_exists) {
#     first_ascending_row <- min(which(profile_data$ascending == 1), na.rm = TRUE)
#
#     if (base_exists) {
#       # Scenario with base and ascending points
#       last_base_row <- max(which(profile_data$base == 1), na.rm = TRUE)
#
#       if (last_base_row < first_ascending_row - 1) {
#         # Mark intermediate points
#         intermediate_rows <- (last_base_row + 1):(first_ascending_row - 1)
#         profile_data$intermediate[intermediate_rows] <- 1
#       }
#     } else {
#       # Scenario 1: No base points, unlabeled points exist before ascending
#       preceding_points <- 1:(first_ascending_row - 1)
#
#       if (length(preceding_points) > 1) {
#         melatonin_values <- profile_data$melatonin[preceding_points]
#         first_ascending_value <- profile_data$melatonin[first_ascending_row]
#
#         # Check threshold condition
#         if (all(melatonin_values < first_ascending_value & melatonin_values < threshold)) {
#           profile_data$intermediate <- 0
#           profile_data$intermediate[preceding_points[-1]] <- 1
#           profile_data$base[preceding_points[1]] <- 1
#         }
#       }
#       # Scenario 2: All ascending except the first point
#       if (first_ascending_row == 2 && profile_data$base[1]==0) {
#         print("here!")
#         # Create a new column 'intermediate', defaulting to 0
#         profile_data <- profile_data %>%
#         dplyr::mutate(intermediate = 0)
#         # Set first point to intermediate
#         profile_data$intermediate[1] <- 1
#         # profile_data$base[1] <- 0
#         # profile_data$slope[1] <- 0
#
#         # Create a new base point 30 minutes earlier
#         new_row <- tibble::tibble(
#           datetime = as.POSIXct(profile_data$datetime[1] - lubridate::minutes(30)),
#           melatonin = profile_data$melatonin[1] / 2,
#           slope = NA,
#           base = 1,
#           ascending = 0,
#           intermediate = 0
#         )
#
#         print(new_row)
#         # Insert the new row at the top
#         profile_data <- dplyr::bind_rows(new_row, profile_data)
#
#         # Remove the 'time' column if it exists
#         if ("time" %in% colnames(profile_data)) {
#           profile_data <- profile_data %>% dplyr::select(-time)
#         }
#       }
#     }
#   }
#
#   # Return the updated profile data
#   return(profile_data)
# }
