# define_ascending_segment <- function(profile_data, threshold = 2.3, interval_limit = lubridate::hours(2)) {
#   # Ensure datetime is sorted
#   profile_data <- profile_data %>% dplyr::arrange(.data$datetime)
#
#   # Identify segments crossing the threshold
#   profile_data <- profile_data %>%
#     dplyr::mutate(
#       transition_to_above = .data$melatonin > threshold &
#         dplyr::lag(.data$melatonin <= threshold, default = FALSE)
#     )
#
#   # Assign groups for each rise above the threshold
#   profile_data <- profile_data %>%
#     dplyr::mutate(
#       rise_group = dplyr::if_else(.data$transition_to_above,
#                                   cumsum(.data$transition_to_above),
#                                   NA_integer_)
#     )
#
#   # Fill rise_group downwards
#   profile_data <- profile_data %>%
#     tidyr::fill(.data$rise_group, .direction = "down")
#
#   # Filter valid rise groups by time intervals
#   rise_times <- profile_data %>%
#     dplyr::filter(!is.na(.data$rise_group)) %>%
#     dplyr::group_by(.data$rise_group) %>%
#     dplyr::summarize(start_time = min(.data$datetime), .groups = "drop")
#
#   rise_times <- rise_times %>%
#     dplyr::mutate(interval = .data$start_time - dplyr::lag(.data$start_time))
#
#   valid_rise_groups <- rise_times %>%
#     dplyr::filter(is.na(.data$interval) | .data$interval >= interval_limit) %>%
#     dplyr::pull(.data$rise_group)
#
#   # Add the ascending column based on valid rise groups
#   profile_data <- profile_data %>%
#     dplyr::mutate(
#       ascending = dplyr::if_else(.data$rise_group %in% valid_rise_groups & .data$melatonin > threshold, 1, 0, 0)
#     )
#
#   # Drop intermediate columns
#   profile_data <- profile_data %>%
#     dplyr::select(-transition_to_above, -rise_group)
#
#   # Return the profile tibble with the ascending column added
#   return(profile_data)
# }

define_ascending_segment <- function(profile_data, threshold = 2.3, interval_limit = lubridate::hours(2)) {
  # Ensure datetime is sorted
  profile_data <- profile_data %>% dplyr::arrange(.data$datetime)

  # Identify segments crossing the threshold
  profile_data <- profile_data %>%
    dplyr::mutate(
      transition_to_above = .data$melatonin > threshold &
        dplyr::lag(.data$melatonin <= threshold, default = FALSE)
    )

  # Assign groups for each rise above the threshold
  profile_data <- profile_data %>%
    dplyr::mutate(
      rise_group = dplyr::if_else(.data$transition_to_above,
                                  cumsum(.data$transition_to_above),
                                  NA_integer_)
    )

  # Fill rise_group downwards
  profile_data <- profile_data %>%
    tidyr::fill(.data$rise_group, .direction = "down")

  # Filter valid rise groups by time intervals
  rise_times <- profile_data %>%
    dplyr::filter(!is.na(.data$rise_group)) %>%
    dplyr::group_by(.data$rise_group) %>%
    dplyr::summarize(start_time = min(.data$datetime), .groups = "drop")

  rise_times <- rise_times %>%
    dplyr::mutate(interval = .data$start_time - dplyr::lag(.data$start_time))

  valid_rise_groups <- rise_times %>%
    dplyr::filter(is.na(.data$interval) | .data$interval >= interval_limit) %>%
    dplyr::pull(.data$rise_group)

  # Add the ascending column based on valid rise groups
  profile_data <- profile_data %>%
    dplyr::mutate(
      ascending = dplyr::if_else(.data$rise_group %in% valid_rise_groups & .data$melatonin > threshold, 1, 0)
    )

  # Add the extra rule: check if the previous point before the first ascending is steep enough
  first_ascending_row <- profile_data %>%
    dplyr::filter(.data$ascending == 1) %>%
    dplyr::slice(1)  # Get the first ascending point

  if (nrow(first_ascending_row) > 0) {
    first_ascending_index <- which(profile_data$datetime == first_ascending_row$datetime)

    # Check if there is a preceding point
    if (first_ascending_index > 1) {
      preceding_row <- profile_data[first_ascending_index - 1, ]
      first_ascending_slope <- first_ascending_row$slope
      preceding_slope <- preceding_row$slope

      # Check the slope rule
      if (abs(preceding_slope) >= abs(first_ascending_slope) / 2) {
        # Mark the previous point as ascending
        profile_data <- profile_data %>%
          dplyr::mutate(
            ascending = dplyr::if_else(.data$datetime == preceding_row$datetime, 1, .data$ascending)
          )
      }
    }
  }

  # Drop intermediate columns
  profile_data <- profile_data %>%
    dplyr::select(-.data$transition_to_above, -.data$rise_group)

  # Return the profile tibble with the ascending column added
  return(profile_data)
}

