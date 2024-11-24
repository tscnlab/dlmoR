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

truncate_ascending_segment <- function(profile_data) {
  # Identify rows that belong to the ascending segment
  ascending_data <- profile_data %>% dplyr::filter(.data$ascending == 1)

  # If no ascending segment exists, issue a warning and return the data unchanged
  if (nrow(ascending_data) == 0) {
    warning("No ascending segments found in the profile data.")
    return(profile_data)  # Return the data unchanged if no ascending segment exists
  }

  # Find the steepest slope in the ascending segment (ignoring NA values)
  max_slope <- max(ascending_data$slope, na.rm = TRUE)

  # Find the index of the last ascending point in the profile data
  last_ascending_index <- max(which(profile_data$ascending == 1))

  # Create a condition for truncating points from the ascending segment
  # Rule 1: Exclude points with a slope <= 0
  # Rule 2: Exclude points with a slope < 50% of the steepest slope (max_slope)
  truncate_condition <- profile_data$ascending == 1 & (
    profile_data$slope <= 0 |  # Rule 1: Non-positive slope
      profile_data$slope < 0.5 * max_slope  # Rule 2: Less than half of the max slope
  )

  # Apply the truncation condition: Set ascending = 0 for points that don't meet the criteria
  profile_data$ascending[truncate_condition] <- 0

  # Explicitly set the last ascending point to 0, regardless of whether it meets the criteria
  profile_data$ascending[last_ascending_index] <- 0

  # Return the modified profile tibble with the truncated ascending segment
  return(profile_data)
}
