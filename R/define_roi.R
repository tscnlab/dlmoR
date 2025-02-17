#' Define Region of Interest (ROI) for Melatonin Profile
#'
#' This function defines the Region of Interest (ROI) in a melatonin profile. The ROI is determined
#' based on the base, intermediate, and ascending segments, as well as the melatonin threshold.
#' The function calculates the horizontal (time) and vertical (melatonin level) bounds of the ROI.
#'
#' @param profile_data A tibble containing the melatonin profile with the following columns:
#'   - `datetime`: A POSIXct column representing timestamps for each measurement.
#'   - `melatonin`: Numeric column representing melatonin concentrations.
#'   - `base`: Binary column (1 for base segment, 0 otherwise).
#'   - `ascending`: Binary column (1 for ascending segment, 0 otherwise).
#'   - `intermediate` (optional): Binary column (1 for intermediate segment, 0 otherwise).
#' @param threshold Numeric. The melatonin concentration threshold used to define the vertical upper bound of the ROI (default = 2.3).
#' @return A list with the following components:
#'   - `x_start`: POSIXct. The start time of the ROI.
#'   - `x_end`: POSIXct. The end time of the ROI.
#'   - `y_min`: Numeric. The lower bound of the ROI (minimum melatonin concentration in the profile).
#'   - `y_max`: Numeric. The upper bound of the ROI (set to the threshold).
#'
#' @examples
#' library(dplyr)
#' library(lubridate)
#'
#' # Example data
#' profile_data <- tibble(
#'   datetime = seq(ymd_hms("2023-01-01 20:00:00"), by = "15 min", length.out = 12),
#'   melatonin = c(1.2, 1.4, 1.5, 2.0, 2.5, 2.3, 2.8, 3.5, 4.2, 4.7, 5.0, 5.2),
#'   base = c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0),
#'   ascending = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1),
#'   intermediate = c(0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0)
#' )
#'
#' # Define the ROI
#' roi <- define_roi(profile_data, threshold = 2.3)
#' print(roi)
#' @export


define_roi <- function(profile_data, threshold = 2.3) {
  # Ensure the profile_data is sorted by datetime for consistent processing
  profile_data <- profile_data %>% dplyr::arrange(datetime)

  # Extract relevant segments
  base_points <- profile_data %>% dplyr::filter(base == 1)
  ascending_points <- profile_data %>% dplyr::filter(ascending == 1)

  # Determine horizontal bounds (x_start and x_end)
  if ("intermediate" %in% colnames(profile_data)) {
    # Case 1: Intermediate segment exists
    intermediate_points <- profile_data %>% dplyr::filter(intermediate == 1)
    x_start <- tail(base_points$datetime, n = 1) +
      0.1 * (intermediate_points$datetime[1] - tail(base_points$datetime, n = 1))
    x_end <- tail(intermediate_points$datetime, n = 1) +
      0.95 * (ascending_points$datetime[1] - tail(intermediate_points$datetime, n = 1))

  } else if (nrow(base_points) < 2) {
    # Case 2: Single-point base segment
    x_start <- base_points$datetime[1] +
      0.05 * (ascending_points$datetime[1] - tail(base_points$datetime, n = 1))
    x_end <- tail(base_points$datetime, n = 1) +
      0.95 * (ascending_points$datetime[1] - tail(base_points$datetime, n = 1))

  } else {
    # Case 3: Multiple base points
    x_start <- decimal_to_posixct(
      mean(posixct_to_decimal(base_points$datetime[(nrow(base_points) - 1):nrow(base_points)], profile_data$datetime)),
      profile_data$datetime
    )
    x_end <- tail(base_points$datetime, n = 1) +
      0.95 * (ascending_points$datetime[1] - tail(base_points$datetime, n = 1))
  }

  # Determine vertical bounds (y_min and y_max)
  y_min <- min(profile_data$melatonin)  # Lowest melatonin concentration in the entire profile
  y_max <- threshold  # Upper bound is set to the threshold

  # Return the ROI as a list
  list(
    x_start = x_start,  # Start time of the ROI
    x_end = x_end,      # End time of the ROI
    y_min = y_min,      # Lower bound of the ROI
    y_max = y_max       # Upper bound of the ROI
  )
}
