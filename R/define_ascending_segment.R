#' Define the ascending segment of a melatonin profile
#'
#' This function identifies and labels the ascending segment of a melatonin profile,
#' where melatonin concentrations rise above a specified threshold and satisfy specific
#' time interval and slope conditions. It also incorporates additional rules to refine
#' the identification of ascending points.
#'
#' @param profile_data A tibble containing melatonin profile data with the following columns:
#'   - `datetime`: A POSIXct column representing timestamps for each measurement.
#'   - `melatonin`: A numeric column representing melatonin concentrations.
#'   - `base`: A binary column (1 for base segments, 0 otherwise) indicating low melatonin baseline segments.
#'   - `slope`: A numeric column representing the rate of change in melatonin concentrations between consecutive points.
#' @param threshold Numeric. The melatonin concentration threshold to define the ascending segment (default = 2.3 pg/mL).
#' @param interval_limit Numeric or lubridate duration. Specifies the minimum time interval between
#' consecutive threshold crossings to consider them as independent ascending events. If numeric, it is interpreted as hours (default = 2 hours).
#' @return A tibble with the following columns:
#'   - `datetime`: The original timestamps.
#'   - `melatonin`: The input melatonin concentrations.
#'   - `base`: The input base column (1 for base segments, 0 otherwise).
#'   - `slope`: The calculated slope of the profile.
#'   - `ascending`: A binary indicator (1 for points in the ascending segment, 0 otherwise).
#'
#' @export

define_ascending_segment <- function(profile_data, threshold = 2.3, interval_limit = lubridate::hours(2)) {
  # Ensure datetime is sorted in ascending order
  profile_data <- profile_data %>% dplyr::arrange(.data$datetime)

  # Identify the melatonin value of the last base point (where base == 1)
  if (any(profile_data$base == 1)) {
    last_base_point <- profile_data %>%
      dplyr::filter(.data$base == 1) %>%
      dplyr::slice_tail(n = 1)

    last_base_melatonin <- last_base_point$melatonin
  } else {
    last_base_point <- dplyr::tibble(datetime = as.POSIXct(NA), melatonin = NA_real_)
    last_base_melatonin <- NA_real_
  }

  # Exclude base segments before threshold checks. If a base segment exists,
  # only consider points after the last base point so early above-threshold
  # samples before the baseline are not treated as the melatonin onset.
  non_base_data <- profile_data %>%
    dplyr::filter(.data$base != 1)

  if (!is.na(last_base_point$datetime)) {
    non_base_data <- non_base_data %>%
      dplyr::filter(.data$datetime > last_base_point$datetime)
  }

  # Identify transitions where melatonin rises above the threshold
  non_base_data <- non_base_data %>%
    dplyr::mutate(
      transition_to_above = .data$melatonin > threshold &
        (dplyr::lag(.data$melatonin, default = last_base_melatonin) <= threshold)
    )

  # Assign groups to each rise above the threshold
  non_base_data <- non_base_data %>%
    dplyr::mutate(
      rise_group = dplyr::if_else(.data$transition_to_above,
                                  cumsum(.data$transition_to_above),
                                  NA_integer_)
    ) %>%
    tidyr::fill(.data$rise_group, .direction = "down") # Propagate group values downwards

  # Calculate the start time of each rise group
  rise_times <- non_base_data %>%
    dplyr::filter(!is.na(.data$rise_group)) %>%
    dplyr::group_by(.data$rise_group) %>%
    dplyr::summarize(start_time = min(.data$datetime), .groups = "drop")

  # Identify valid rise groups based on the interval limit. Within a short
  # interval, repeated threshold crossings are treated as the same rise event,
  # and the later crossing is used to represent that event.
  valid_rise_groups <- integer(0)
  pending_rise_group <- NULL
  previous_start_time <- NULL
  interval_limit_seconds <- as.numeric(interval_limit)

  for (i in seq_len(nrow(rise_times))) {
    current_start_time <- rise_times$start_time[i]
    current_rise_group <- rise_times$rise_group[i]

    if (is.null(pending_rise_group)) {
      pending_rise_group <- current_rise_group
    } else if (as.numeric(difftime(current_start_time, previous_start_time, units = "secs")) <= interval_limit_seconds) {
      pending_rise_group <- current_rise_group
    } else {
      valid_rise_groups <- c(valid_rise_groups, pending_rise_group)
      pending_rise_group <- current_rise_group
    }

    previous_start_time <- current_start_time
  }

  if (!is.null(pending_rise_group)) {
    valid_rise_groups <- c(valid_rise_groups, pending_rise_group)
  }

  if (length(valid_rise_groups) > 0) {
    valid_rise_groups <- valid_rise_groups[1]
  }

  # Mark rows as ascending if they belong to valid rise groups and exceed the threshold
  non_base_data <- non_base_data %>%
    dplyr::mutate(
      ascending = dplyr::if_else(.data$rise_group %in% valid_rise_groups & .data$melatonin > threshold, 1, 0)
    )

  # Merge the updated ascending column back into the original profile_data
  profile_data <- profile_data %>%
    dplyr::left_join(non_base_data %>% dplyr::select(datetime, ascending), by = "datetime") %>%
    dplyr::mutate(ascending = dplyr::coalesce(.data$ascending, 0)) # Fill NA with 0 for non-ascending rows

  # Use accepted ascending points to identify steepest slopes. This avoids
  # re-introducing threshold crossings that were rejected by interval_limit.
  accepted_ascending_data <- profile_data %>%
    dplyr::filter(.data$ascending == 1)

  # Identify the steepest slope among accepted ascending rows
  steepest_slope <- if (all(is.na(accepted_ascending_data$slope))) NA_real_ else max(accepted_ascending_data$slope, na.rm = TRUE)

  # Identify segments with slopes >= half of the steepest slope
  steep_segments <- accepted_ascending_data %>%
    dplyr::filter(.data$slope >= steepest_slope / 2) %>%
    dplyr::pull(.data$datetime)

  # Include segments between the steepest points as part of the ascending segment
  if (length(steep_segments) > 1) {
    min_datetime <- min(steep_segments)
    max_datetime <- max(steep_segments)

    profile_data <- profile_data %>%
      dplyr::mutate(
        ascending = dplyr::if_else(
          (.data$datetime >= min_datetime & .data$datetime <= max_datetime & .data$base != 1) | .data$ascending == 1,
          1,
          0
        )
      )
  }

  # Ensure the row before the first ascending point is included if the slope condition is met
  first_ascending_row <- profile_data %>%
    dplyr::filter(.data$ascending == 1) %>%
    dplyr::slice(1)

  if (nrow(first_ascending_row) > 0) {
    first_ascending_index <- dplyr::first(which(profile_data$datetime == first_ascending_row$datetime))

    if (nrow(last_base_point) > 0) {
      last_base_index <- dplyr::first(which(profile_data$datetime == last_base_point$datetime))
    } else {
      last_base_index <- NA_integer_
    }

    if (first_ascending_index > 1) {
      current_index <- first_ascending_index - 1
      keep_checking <- TRUE
      rows_to_update <- c()

      while (current_index > 0 && keep_checking) {
        preceding_row <- profile_data[current_index, ]
        preceding_row_index <- dplyr::first(which(profile_data$datetime == preceding_row$datetime))
        first_ascending_slope <- first_ascending_row$slope
        preceding_slope <- preceding_row$slope

        if (!is.na(last_base_index) && length(last_base_index) > 0) {
          base_condition <- (preceding_row_index != last_base_index)
        } else {
          base_condition <- TRUE
        }

        if (!is.na(preceding_slope) && (preceding_slope >= first_ascending_slope / 2) && base_condition) {
          rows_to_update <- c(rows_to_update, preceding_row_index)
          current_index <- current_index - 1
        } else {
          keep_checking <- FALSE
        }
      }

      profile_data <- profile_data %>%
        dplyr::mutate(ascending = dplyr::if_else(.data$datetime %in% profile_data$datetime[rows_to_update], 1, .data$ascending))
    }
  }

  profile_data <- profile_data %>%
    dplyr::select(-tidyselect::any_of(c("transition_to_above", "rise_group")))

  return(profile_data)
}
