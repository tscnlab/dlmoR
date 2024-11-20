#' Data Validation and Preprocessing Functions
#'
#' This file contains functions for validating and preprocessing profile data for analysis.
#' These functions help ensure that profiles meet minimum requirements for rise, length,
#' and ending criteria before proceeding with further analysis.
#'
#' ## Functions included:
#' - `check_rise_threshold()`: Checks if the profile has a sustained rise above a specified threshold.
#' - `check_data_points()`: Ensures the profile has a minimum number of data points.
#' - `trim_end_below_threshold()`: Trims nodes below the threshold at the end of the profile.
#' - `preprocess_profile()`: Wrapper function that combines all checks and preprocessing steps.
#'
#' @param profile A numeric vector representing the profile data points.
#' @param threshold Numeric value for the minimum rise threshold, default is 2.3.
#' @param min_increase_points Integer specifying the minimum number of consecutive points
#'   above the threshold to establish a rise, default is 3.
#' @param min_points Integer specifying the minimum number of data points required for the profile,
#'   default is 3.
#' @return Depending on the function, returns either a boolean (for validation functions) or a
#'   numeric vector (for the processed profile). If criteria are not met, an error message is returned.
#'
#' ## Example Usage
#' ```R
#' profile <- c(1.2, 2.1, 2.4, 2.5, 2.8, 2.9, 1.9, 2.2)
#'
#' # Check if profile has a significant rise
#' check_rise_threshold(profile, threshold = 2.3, min_increase_points = 3)
#'
#' # Check if profile has sufficient data points
#' check_data_points(profile, min_points = 3)
#'
#' # Trim nodes below threshold at the end of the profile
#' trimmed_profile <- trim_end_below_threshold(profile, threshold = 2.3)
#'
#' # Full preprocessing, with all checks and trimming
#' processed_profile <- preprocess_profile(profile, threshold = 2.3, min_increase_points = 3, min_points = 3)
#' ```
#'
#' @name validate_melatonin_profile
NULL

# validate_melatonin_profile(profile, threshold = 2.3, min_increase_points = 3){
#
# }


#' @inheritParams validate_melatonin_profile
check_rise_threshold <- function(profile, threshold = 2.3, min_increase_points = 3) {
  library(dplyr)
  library(tibble)

  profile_tbl <- tibble(value = profile) %>%
    mutate(
      above_threshold = value >= threshold,
      # Handle the first row explicitly: prev_below_threshold should be FALSE for the first row
      prev_below_threshold = ifelse(row_number() == 1, FALSE, lag(value < threshold)),
      # Transition occurs only if crossing from below to above threshold
      transition = above_threshold & prev_below_threshold
    ) %>%
    # Group consecutive above-threshold points starting from valid transitions
    mutate(group = cumsum(transition))

  # Print profile_tbl to see intermediate results
  print("Profile Table After Mutation:")
  print(profile_tbl,n=46)  # Debugging: Inspect intermediate results

  # Group by `group` and calculate sustained length
  summarized_tbl <- profile_tbl %>%
    group_by(group) %>%
    summarize(sustained_length = sum(above_threshold), .groups = "drop")

  # Print the grouped summary (sustained length per group)
  print("Summarized Sustained Length Per Group:")
  print(summarized_tbl)

  # Final check for any sustained rise, excluding group 0 (no valid transition)
  final_result <- summarized_tbl %>%
    filter(group != 0) %>%  # Filter out group 0
    summarize(any_rise = any(sustained_length >= min_increase_points))

  # Print the final result for any_rise
  print("Final Result (any_rise):")
  print(final_result)

  return(final_result$any_rise)  # Return the result
}




#' @inheritParams validate_melatonin_profile
check_data_length <- function(profile, min_points = 3) {
  length(profile) >= min_points
  data_length<-length(profile)
  print(data_length)
}

#' @inheritParams validate_melatonin_profile
trim_end_below_threshold <- function(profile, threshold = 2.3) {
  last_above_threshold <- max(which(profile >= threshold))
  profile[1:last_above_threshold]
}

#' @export
#' @inheritParams validate_melatonin_profile
preprocess_profile <- function(profile, threshold = 2.3, min_increase_points = 3, min_points = 3) {
  if (!check_data_length(profile, min_points)) {
    stop("Insufficient number of data points for analysis")
  }
  if (!check_rise_threshold(profile, threshold, min_increase_points)) {
    stop("No significant rise detected in profile for analysis")
  }
  trimmed_profile <- trim_end_below_threshold(profile, threshold)

  return(trimmed_profile)
}
