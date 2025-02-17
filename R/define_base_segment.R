#' Define the base segment of a melatonin profile
#'
#' This function identifies and labels the base segment of a melatonin profile, where melatonin
#' concentrations are below a specified threshold and the slope of the profile is non-positive
#' (indicating stable or decreasing melatonin levels).
#'
#' @param profile_tibble A tibble containing `datetime`, `melatonin`, and `time` columns.
#'   - `datetime`: A POSIXct column representing timestamps for each melatonin measurement.
#'   - `melatonin`: A numeric column representing melatonin concentrations at each time point.
#'   - `time`: A numeric column representing the time in hh:mm:ss format
#' @param threshold Numeric. The melatonin concentration threshold to define the base segment (default = 2.3 pg/mL).
#' @return A tibble with the base segment labeled (1 for base, 0 for non-base), including:
#'   - `datetime`: The original timestamps.
#'   - `melatonin`: The input melatonin concentrations.
#'   - `time`: The hh:mm:ss time format
#'   - `slope`: The rate of change in melatonin concentrations between consecutive time points, calculated as
#'     \deqn{(melatonin[i+1] - melatonin[i]) / (time[i+1] - time[i])}.
#'   - `base`: A binary indicator (1 for points in the base segment, 0 otherwise).
#'
#' The slope helps determine periods when melatonin levels are stable or decreasing, which are key characteristics
#' of the base segment.
#'
#' @examples
#' # Example data
#' library(dplyr)
#' library(lubridate)
#'
#' profile_data <- tibble(
#'   datetime = seq(ymd_hms("2023-01-01 20:00:00"), by = "15 min", length.out = 12),
#'   melatonin = c(1.2, 1.4, 1.5, 1.7, 2.0, 2.3, 2.8, 3.5, 4.2, 4.7, 5.0, 5.2),
#'   time = decimal_date(seq(ymd_hms("2023-01-01 20:00:00"), by = "15 min", length.out = 12))
#' )
#'
#' # Define the base segment using a threshold of 2.3
#' base_segment <- define_base_segment(profile_tibble = profile_data, threshold = 2.3)
#'
#' # View the results
#' print(base_segment)
#'
#' @export


define_base_segment <- function(profile_tibble, threshold = 2.3) {
  # Extract melatonin values from the input tibble (profile_tibble)
  profile_tibble <- profile_tibble %>% dplyr::arrange(.data$datetime)
  profile <- profile_tibble$melatonin
  decimal_time <- posixct_to_decimal(profile_tibble$datetime, profile_tibble$datetime)
  # Calculate slopes based on the melatonin profile
  slopes <- calculate_slopes(profile, decimal_time)

  # Create a tibble for the full profile with slopes and base labels
  profile_data <- profile_tibble %>%
    dplyr::mutate(
      slope = c(NA, slopes),                     # Append NA at start to align slope with profile length; assigns slope value to the righthand node of a segment
      base = rep(0, length(profile))              # Initialize base column to 0 for all points
    )

  # Identify base segment based on slope and threshold
  profile_data <- profile_data %>%
    dplyr::mutate(
      base = dplyr::if_else(
        .data$slope <= 0 & (.data$melatonin <= threshold | dplyr::lead(.data$melatonin, default = NA) <= threshold),# second part of OR statement, is to identify downward threshold crossing
        1,  # Mark base as 1
        .data$base,  # Keep current base value
        missing = .data$base  # Keep NA for missing values
      )
    ) %>%
    # Propagate the base segment back from the rightmost identified base
    dplyr::mutate(
      base = dplyr::if_else(dplyr::row_number() <= max(which(.data$base == 1), na.rm = TRUE) , 1, .data$base, missing = .data$base)
    )

  # Return the full tibble with datetime, melatonin, time, slope, and base columns
  return(profile_data)
}
