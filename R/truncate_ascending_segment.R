# truncate_ascending_segment <- function(profile_data) {
#   # Identify rows that belong to the ascending segment
#   ascending_data <- profile_data %>% dplyr::filter(.data$ascending == 1)
#
#   # If no ascending segment exists, issue a warning and return the data unchanged
#   if (nrow(ascending_data) == 0) {
#     warning("No ascending segments found in the profile data.")
#     return(profile_data)  # Return the data unchanged if no ascending segment exists
#   }
#
#   # Find the steepest slope in the ascending segment (ignoring NA values)
#   max_slope <- max(ascending_data$slope, na.rm = TRUE)
#
#   # Identify the last ascending index
#   last_ascending_index <- max(which(profile_data$ascending == 1))
#
#   # Check the validity of the current point based on rules and exclude from the NEXT point onward if invalid
#   profile_data <- profile_data %>%
#     dplyr::mutate(
#       # Rule 1: slope > 0
#       # Rule 2: slope >= 50% of max_slope
#       valid_ascending = .data$ascending == 1 &
#         (.data$slope > 0 & .data$slope >= 0.5 * max_slope)
#     )
#
#   # Identify invalid indices and exclude from the NEXT point onward
#   invalid_indices <- which(!profile_data$valid_ascending & profile_data$ascending == 1)
#   if (length(invalid_indices) > 0) {
#     # Find the first invalid point
#     first_invalid_index <- min(invalid_indices)
#
#     # Set all points AFTER the invalid point to not ascending
#     profile_data <- profile_data %>%
#       dplyr::mutate(
#         ascending = dplyr::if_else(
#           dplyr::row_number() > first_invalid_index,
#           0,
#           .data$ascending
#         )
#       )
#   }
#
#   # Drop the temporary valid_ascending column
#   profile_data <- profile_data %>% dplyr::select(-valid_ascending)
#
#   # Return the modified profile tibble with the truncated ascending segment
#   return(profile_data)
# }
#
truncate_ascending_segment <- function(profile_data) {
  # Ensure there is an ascending segment to work with
  ascending_data <- profile_data %>% dplyr::filter(.data$ascending == 1)

  if (nrow(ascending_data) == 0) {
    warning("No ascending segments found in the profile data.")
    return(profile_data)  # Return unchanged if no ascending segment
  }

  # Helper function to check if rules are satisfied
  check_rules <- function(profile_data, max_slope) {
    # Get the ascending segment
    ascending_data <- profile_data %>% dplyr::filter(.data$ascending == 1)

    if (nrow(ascending_data) < 2) {
      return(TRUE)  # If no ascending segment or single point, rules are satisfied
    }

    # Extract the two rightmost points in the ascending segment
    rightmost_points <- ascending_data %>%
      dplyr::arrange(desc(dplyr::row_number())) %>%
      head(2)

    # Calculate the slope of the rightmost segment
    rightmost_slope <- (rightmost_points$melatonin[2] - rightmost_points$melatonin[1]) /
      (as.numeric(difftime(rightmost_points$datetime[2], rightmost_points$datetime[1], units = "secs")))

    # Rule (1a) and (1b): Rightmost slope cannot be zero or negative
    if (rightmost_slope <= 0) {
      return(FALSE)
    }

    # Rule (2a) and (2b): Parallelogram diagonal slope rules
    diag1_slope <- (rightmost_points$melatonin[2] - ascending_data$melatonin[1]) /
      (as.numeric(difftime(rightmost_points$datetime[2], ascending_data$datetime[1], units = "secs")))
    diag2_slope <- (ascending_data$melatonin[nrow(ascending_data)] - ascending_data$melatonin[1]) /
      (as.numeric(difftime(ascending_data$datetime[nrow(ascending_data)], ascending_data$datetime[1], units = "secs")))

    diag_ratio <- abs(diag1_slope / diag2_slope)
    if (diag_ratio < 0.5 || diag_ratio < 0) {
      return(FALSE)
    }

    # Rule (3): Rightmost slope cannot be < 1/2 of the steepest ascending slope
    if (rightmost_slope < 0.5 * max_slope) {
      return(FALSE)
    }

    return(TRUE)
  }

  # Calculate the maximum slope in the ascending segment
  max_slope <- max(profile_data$slope[profile_data$ascending == 1], na.rm = TRUE)

  # Iteratively truncate until all rules are satisfied
  while (!check_rules(profile_data, max_slope)) {
    # Identify the last ascending index
    last_ascending_index <- max(which(profile_data$ascending == 1))

    # Set the last point in the ascending segment to not ascending
    profile_data <- profile_data %>%
      dplyr::mutate(
        ascending = dplyr::if_else(
          dplyr::row_number() == last_ascending_index,
          0,
          .data$ascending
        )
      )
  }

  return(profile_data)
}
