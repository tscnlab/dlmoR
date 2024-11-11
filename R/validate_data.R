# Check if there is a meaningful rise >= threshold
#' Title
#'
#' @param profile
#' @param threshold
#' @param min_increase_points
#'
#' @return
#' @export
#'
#' @examples
check_rise_threshold <- function(profile, threshold = 2.3, min_increase_points = 3) {
  # Find indices where the profile crosses the threshold
  above_threshold <- which(profile >= threshold)

  # If no points are above the threshold, return an error
  if (length(above_threshold) == 0) {
    stop("Error: Profile does not have a significant rise above the threshold.")
  }

  # Check for consecutive increases leading up to a threshold-crossing point
  for (i in 2:length(profile)) {
    if (profile[i] >= threshold && sum(diff(profile[1:i]) > 0) >= min_increase_points) {
      return(TRUE)
    }
  }

  # If no meaningful rise is found
  stop("Error: Profile lacks a sustained rise above the threshold, indicating no dynamic part.")
}

# Check for sufficient number of data points
check_data_points <- function(profile, min_points = 3) {
  if (length(profile) < min_points) {
    stop("Error: Insufficient data points for analysis (minimum required is 3).")
  }
  return(TRUE)
}

# Trim nodes below threshold at the end of the profile
#' Title
#'
#' @param profile
#' @param threshold
#'
#' @return
#' @export
#'
#' @examples
trim_end_below_threshold <- function(profile, threshold = 2.3) {
  last_above_threshold <- max(which(profile >= threshold), na.rm = TRUE)

  # If trimming occurs, return a warning
  if (last_above_threshold < length(profile)) {
    warning("Warning: Trimming nodes below threshold at the end of the profile.")
  }

  return(profile[1:last_above_threshold])
}

# Main function
#' Title
#'
#' @param data
#' @param threshold
#' @param min_increase_points
#' @param min_points
#'
#' @return
#' @export
#'
#' @examples
preprocess_profiles <- function(data, threshold = 2.3, min_increase_points = 3, min_points = 3) {
  processed_data <- list()

  for (i in seq_along(data)) {
    profile <- data[[i]]

    # Check for meaningful rise and handle errors
    check_rise_threshold(profile, threshold, min_increase_points)

    # Check for sufficient data points and handle errors
    check_data_points(profile, min_points)

    # Trim end nodes below threshold
    trimmed_profile <- trim_end_below_threshold(profile, threshold)

    # Append to processed data
    processed_data[[length(processed_data) + 1]] <- trimmed_profile
  }

  return(processed_data)
}

# Example usage
# data <- list(c(1.2, 1.9, 2.1, 2.5, 3.0), c(2.1, 1.9), c(1.8, 2.0, 2.4, 2.7, 2.6, 2.2))
# processed_data <- preprocess_profiles(data, threshold = 2.3, min_increase_points = 3, min_points = 3)
# print(processed_data)
