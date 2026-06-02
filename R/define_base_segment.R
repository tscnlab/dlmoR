#' Define the base segment of a melatonin profile
#'
#' This function identifies and labels the base segment of a melatonin profile, where melatonin
#' concentrations are below a specified threshold and the slope of the profile is non-positive
#' (indicating stable or decreasing melatonin levels).
#'
#' @param profile_tibble A tibble containing `datetime`, `melatonin`, and `time` columns.
#'   - `datetime`: A POSIXct column representing timestamps for each melatonin measurement.
#'   - `melatonin`: A numeric column representing melatonin concentrations at each time point.
#'   - `time`: A numeric column representing the time in hh:mm:ss format.
#' @param threshold Numeric. The melatonin concentration threshold to define the base segment (default = 2.3 pg/mL).
#' @return A tibble with the base segment labeled (1 for base, 0 for non-base), including:
#'   - `datetime`: The original timestamps.
#'   - `melatonin`: The input melatonin concentrations.
#'   - `time`: The hh:mm:ss time format.
#'   - `slope`: The rate of change in melatonin concentrations between consecutive time points.
#'   - `base`: A binary indicator (1 for points in the base segment, 0 otherwise).
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
  # Ensure the data is sorted by datetime
  profile_tibble <- profile_tibble %>% dplyr::arrange(.data$datetime)

  # Extract melatonin values
  profile <- profile_tibble$melatonin
  decimal_time <- posixct_to_decimal(profile_tibble$datetime, profile_tibble$datetime)

  # Calculate slopes based on melatonin profile
  slopes <- calculate_slopes(profile, decimal_time)

  # Create a tibble with slopes and base labels
  profile_data <- profile_tibble %>%
    dplyr::mutate(
      slope = c(NA, slopes),  # Append NA at start to align slope with profile length
      base = rep(0, length(profile))  # Initialize base column to 0 for all points
    )

  # Identify base segment based on slope and threshold
  profile_data <- profile_data %>%
    dplyr::mutate(
      base = dplyr::if_else(
        .data$slope <= 0 & (.data$melatonin <= threshold | dplyr::lead(.data$melatonin, default = NA) <= threshold),
        1,  # Mark as base
        .data$base,  # Keep existing values
        missing = .data$base  # Keep NA values
      )
    )

  # **Fix for max() issue: Handle cases where there are no base points**
  last_base_index <- if (any(profile_data$base == 1)) {
    max(which(profile_data$base == 1), na.rm = TRUE)
  } else {
    NA_integer_  # Ensure NA instead of -Inf
  }

  # Propagate base values correctly
  profile_data <- profile_data %>%
    dplyr::mutate(
      base = dplyr::if_else(
        dplyr::row_number() <= last_base_index & .data$melatonin <= threshold,
        1,
        .data$base,
        missing = .data$base
      )
    )

  # Keep the base segment anchored to the initial low-melatonin portion of
  # the profile. A later dip below threshold after a sustained rise should not
  # move the end of the base segment forward. Single above-threshold blips are
  # not treated as the start of the rise.
  first_rise_index <- profile_data %>%
    dplyr::mutate(
      row_index = dplyr::row_number(),
      above_threshold = .data$melatonin > threshold,
      transition_to_above = .data$above_threshold &
        dplyr::lag(.data$melatonin <= threshold, default = FALSE),
      rise_group = dplyr::if_else(.data$transition_to_above, cumsum(.data$transition_to_above), NA_integer_)
    ) %>%
    tidyr::fill(.data$rise_group, .direction = "down") %>%
    dplyr::filter(.data$above_threshold, !is.na(.data$rise_group)) %>%
    dplyr::group_by(.data$rise_group) %>%
    dplyr::summarize(
      row_index = dplyr::first(.data$row_index),
      n_above_threshold = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::filter(.data$n_above_threshold >= 2) %>%
    dplyr::slice(1) %>%
    dplyr::pull(.data$row_index)

  if (length(first_rise_index) > 0) {
    profile_data <- profile_data %>%
      dplyr::mutate(
        base = dplyr::if_else(dplyr::row_number() >= first_rise_index, 0, .data$base)
      )
  }

  # Return updated tibble
  return(profile_data)
}
