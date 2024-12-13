# truncate_ascending_segment <- function(profile_data) {
#   # Identify rows that belong to the ascending segment
#   ascending_data <- profile_data %>% dplyr::filter(.data$ascending == 1)
#
#   if (nrow(ascending_data) == 0) {
#     warning("No ascending segments found in the profile data.")
#     return(profile_data)  # If there are no ascending segments, return the data as it is
#   }
#
#   # Find the steepest slope in the ascending segment (ignoring NA values)
#   max_slope <- max(ascending_data$slope, na.rm = TRUE)
#
#   # Find the index of the last ascending point
#   last_ascending_index <- max(which(profile_data$ascending == 1))
#
#   # Check slopes from the second to last ascending segment
#   # We create a condition to truncate points that don't meet the criteria
#   truncate_condition <- profile_data$ascending == 1 & (
#     profile_data$slope <= 0 | profile_data$slope < 0.5 * max_slope
#   )
#
#   # Apply the condition to set ascending = 0 for points that should be truncated
#   profile_data$ascending[truncate_condition] <- 0
#
#   # Explicitly set the last ascending point to 0, even if it meets the criteria
#   profile_data$ascending[last_ascending_index] <- 0
#
#   # Return the modified profile tibble with the truncated ascending segment
#   return(profile_data)
# }

# last working version (ish)
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
#   # Find the index of the last ascending point in the profile data
#   last_ascending_index <- max(which(profile_data$ascending == 1))
#
#   # Create a condition for truncating points from the ascending segment
#   # Rule 1: Exclude points with a slope <= 0
#   # Rule 2: Exclude points with a slope < 50% of the steepest slope (max_slope)
#   truncate_condition <- profile_data$ascending == 1 & (
#     profile_data$slope <= 0 |  # Rule 1: Non-positive slope
#       profile_data$slope < 0.5 * max_slope  # Rule 2: Less than half of the max slope
#   )
#
#   # Apply the truncation condition: Set ascending = 0 for points that don't meet the criteria
#   profile_data$ascending[truncate_condition] <- 0
#
#   # Explicitly set the last ascending point to 0, regardless of whether it meets the criteria
#   profile_data$ascending[last_ascending_index] <- 0
#
#   # Return the modified profile tibble with the truncated ascending segment
#   return(profile_data)
# }

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
#   # Create a condition to check the truncation rules for all ascending points
#   # Rule 1: Slope must be > 0
#   # Rule 2: Slope must be >= 50% of the max slope
#   valid_ascending <- profile_data$slope > 0 & profile_data$slope >= 0.40 * max_slope
#
#   # Ensure only the rightmost ascending segments are evaluated for truncation
#   # Create a mask for indices after the last valid ascending point
#   valid_indices <- which(valid_ascending & profile_data$ascending == 1)
#   if (length(valid_indices) > 0) {
#     # Determine the last valid index based on truncation rules
#     last_valid_index <- max(valid_indices)
#
#     # Set all points after the last valid point to not ascending
#     profile_data$ascending[(last_valid_index + 1):last_ascending_index] <- 0
#   }
#
#   # Explicitly set the last ascending point to 0
#   profile_data$ascending[last_ascending_index] <- 0
#
#   # Return the modified profile tibble with the truncated ascending segment
#   return(profile_data)
# }

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
#   # Identify valid ascending points based on the truncation rules:
#   # Rule 1: Slope must be > 0
#   # Rule 2: Slope must be >= 50% of the max slope
#   profile_data <- profile_data %>%
#     dplyr::mutate(
#       valid_ascending = .data$slope > 0 & .data$slope >= 0.5 * max_slope & .data$ascending == 1
#     )
#
#   # Find the last valid ascending point that satisfies the conditions
#   last_valid_index <- max(which(profile_data$valid_ascending), na.rm = TRUE)
#
#   # Set all points AFTER the last valid ascending point to 0
#   profile_data <- profile_data %>%
#     dplyr::mutate(
#       ascending = dplyr::if_else(
#         dplyr::row_number() > last_valid_index & dplyr::row_number() <= last_ascending_index,
#         0,
#         .data$ascending
#       )
#     )
#
#   # Explicitly set the last ascending point to 0
#   profile_data <- profile_data %>%
#     dplyr::mutate(
#       ascending = dplyr::if_else(
#         dplyr::row_number() == last_ascending_index,
#         0,
#         .data$ascending
#       )
#     )
#
#   # Drop the temporary valid_ascending column
#   profile_data <- profile_data %>%
#     dplyr::select(-valid_ascending)
#
#   # Return the modified profile tibble with the truncated ascending segment
#   return(profile_data)
# }
#

truncate_ascending_segment_orig <- function(profile_data) { #TODO THIS SCRIPT CURRENTLY DOES NOT TRUNCATE PARALLELOGRAMMMM; merge with wip script
 print("old")
   # Identify rows that belong to the ascending segment
  ascending_data <- profile_data %>% dplyr::filter(.data$ascending == 1)

  # If no ascending segment exists, issue a warning and return the data unchanged
  if (nrow(ascending_data) == 0) {
    warning("No ascending segments found in the profile data.")
    return(profile_data)  # Return the data unchanged if no ascending segment exists
  }

  # Find the steepest slope in the ascending segment (ignoring NA values)
  max_slope <- max(ascending_data$slope, na.rm = TRUE)

  # Identify the last ascending index
  last_ascending_index <- max(which(profile_data$ascending == 1))

  # Check the validity of the current point based on rules and exclude from the NEXT point onward if invalid
  profile_data <- profile_data %>%
    dplyr::mutate(
      # Rule 1: slope > 0
      # Rule 2: slope >= 50% of max_slope
      valid_ascending = .data$ascending == 1 &
        (.data$slope > 0 & .data$slope >= 0.5 * max_slope)
    )

  # Identify invalid indices and exclude from the NEXT point onward
  invalid_indices <- which(!profile_data$valid_ascending & profile_data$ascending == 1)
  if (length(invalid_indices) > 0) {
    # Find the first invalid point
    first_invalid_index <- min(invalid_indices)

    # Set all points AFTER the invalid point to not ascending
    profile_data <- profile_data %>%
      dplyr::mutate(
        ascending = dplyr::if_else(
          dplyr::row_number() > first_invalid_index,
          0,
          .data$ascending
        )
      )
  }

  # Drop the temporary valid_ascending column
  profile_data <- profile_data %>% dplyr::select(-valid_ascending)

  # Return the modified profile tibble with the truncated ascending segment
  #print("post-truncation")
  #print(profile_data, n = 28)
  return(profile_data)
}

