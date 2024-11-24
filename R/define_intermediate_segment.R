define_intermediate_segment <- function(profile_data) {
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
