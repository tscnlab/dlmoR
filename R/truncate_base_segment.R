truncate_base_segment <- function(profile_data, threshold = 2.3) {
  # Ensure the profile is sorted by datetime
  profile_data <- profile_data %>% dplyr::arrange(.data$datetime)

  # Truncate points in the base segment that are above the threshold
  profile_data <- profile_data %>%
    dplyr::mutate(
      base = dplyr::if_else(.data$base == 1 & .data$melatonin > threshold, 0, .data$base)
    )

  # Check if all remaining segments belong to the ascending part
  if (all(profile_data$base == 0)) {
    warning("Hockey-stick time: no base part")
    return(profile_data)  # Return the data as is with base set to 0
  }

  # Return the modified profile tibble
  return(profile_data)
}
