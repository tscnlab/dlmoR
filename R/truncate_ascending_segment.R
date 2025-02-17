#' Truncate Ascending Segment of Melatonin Profile
#'
#' This function truncates the ascending segment of a melatonin profile based on a set of
#' rules to ensure valid slopes and segment consistency. It iteratively modifies the
#' `ascending` column of the input data to satisfy the rules and uses the
#' `parallelogram_fit` function as a diagnostic measure to validate the final segment. The
#' `parallelogram_fit` function evaluates the geometric structure of the ascending
#' segment and ensures the segment satisfies predefined geometric constraints, including
#' slope ratios of lateral and diagonal shape segments. This ensures that only the steadily
#' increasing part of the melatonin rise is taken into account when determining the DLMO point
#' and that slower rises, or drops in melatonin levels are not.
#'
#'
#' @param profile_data A tibble containing the melatonin profile with the following columns:
#'   - `datetime`: A POSIXct column with timestamps.
#'   - `melatonin`: Numeric column representing melatonin concentrations.
#'   - `slope`: Numeric column with the slope between consecutive melatonin points.
#'   - `ascending`: Binary column (1 for ascending segment, 0 otherwise).
#' @return A list containing:
#'   - `profile`: The updated profile_data tibble with truncated ascending segments.
#'   - `plll`: The result of the parallelogram fit diagnostic (from `parallelogram_fit`).
#'
#'
#' This function ensures the following rules are satisfied:
#' 1. The rightmost slope in the ascending segment cannot be less than half the steepest slope in the segment.
#' 2. The ascending segment must adhere to additional conditions validated by the `parallelogram_fit` function.
#'
#' @examples
#' library(dplyr)
#' library(lubridate)
#'
#' # Example data
#' profile_data <- tibble(
#'   datetime = seq(ymd_hms("2023-01-01 20:00:00"), by = "15 min", length.out = 12),
#'   melatonin = c(1.2, 1.4, 1.5, 1.7, 2.0, 2.3, 2.8, 3.5, 4.2, 4.7, 5.0, 5.2),
#'   slope = c(NA, diff(c(1.2, 1.4, 1.5, 1.7, 2.0, 2.3, 2.8, 3.5, 4.2, 4.7, 5.0, 5.2))),
#'   ascending = c(0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1)
#' )
#'
#' # Truncate the ascending segment
#' result <- truncate_ascending_segment(profile_data)
#' print(result$profile)
#' print(result$plll)
#' @export


truncate_ascending_segment <- function(profile_data) {
  # Ensure there is an ascending segment to work with
  ascending_data <- profile_data %>% dplyr::filter(.data$ascending == 1)

  # If no ascending segment exists, issue a warning and return the data unchanged
  if (nrow(ascending_data) == 0) {
    warning("No ascending segments found in the profile data.")
    return(profile_data)
  }

  # Find the steepest slope in the ascending segment
  max_slope <- max(profile_data$slope[profile_data$ascending == 1], na.rm = TRUE)

  # Rule Validation Loop: Iteratively truncate until all rules are satisfied
  while (!check_rules(profile_data, max_slope)) {
    # Identify the last ascending index
    last_ascending_index <- max(which(profile_data$ascending == 1))

    # Remove the last point in the ascending segment
    profile_data <- profile_data %>%
      dplyr::mutate(
        ascending = dplyr::if_else(
          dplyr::row_number() == last_ascending_index,
          0,
          .data$ascending
        )
      )
  }

  # Rule (3) Validation: Ensure the rightmost slope is valid
  while (!check_last(profile_data)) {
    # Identify the last ascending index
    last_ascending_index <- max(which(profile_data$ascending == 1))

    # Remove the last point in the ascending segment
    profile_data <- profile_data %>%
      dplyr::mutate(
        ascending = dplyr::if_else(
          dplyr::row_number() == last_ascending_index,
          0,
          .data$ascending
        )
      )
  }

  # Perform a parallelogram fit diagnostic on the truncated profile
  plll <- parallelogram_fit(profile_data)

  # Return the updated profile and parallelogram fit diagnostic
  return(list(profile = profile_data, plll = plll))
}

# Helper Function: Check if the rightmost slope is valid
check_last <- function(profile_data) {
  # Find the steepest slope in the ascending segment
  max_slope <- max(profile_data$slope[profile_data$ascending == 1], na.rm = TRUE)

  # Identify the last ascending index
  last_ascending_index <- max(which(profile_data$ascending == 1))

  # Check if the rightmost slope satisfies the rule
  return(profile_data$slope[last_ascending_index] >= 0.5 * max_slope)
}

# Helper Function: Check if all rules are satisfied
check_rules <- function(profile_data, max_slope) {
  # Extract the ascending segment
  ascending_data <- profile_data %>% dplyr::filter(.data$ascending == 1)

  # Rule: If fewer than two points, rules are satisfied
  if (nrow(ascending_data) < 2) {
    return(TRUE)
  }

  # Calculate the slope of the rightmost segment
  rightmost_points <- ascending_data %>%
    dplyr::arrange(desc(dplyr::row_number())) %>%
    head(2)
  rightmost_slope <- (rightmost_points$melatonin[2] - rightmost_points$melatonin[1]) /
    (as.numeric(difftime(rightmost_points$datetime[2], rightmost_points$datetime[1], units = "secs")))

  # Rule (1): Rightmost slope must be positive
  if (rightmost_slope <= 0) {
    return(FALSE)
  }

  # Rule (2): Parallelogram fit must be valid
  plll <- parallelogram_fit(profile_data)
  if (plll$flag) {
    return(FALSE)
  }

  # If all rules are satisfied
  return(TRUE)
}

