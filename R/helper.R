# calculate slopes between data points
calculate_slopes <- function(profile, decimal_time) {
  # Ensure datetime is sorted
  # Ensure datetime is sorted
  # profile_data <- profile_data %>% dplyr::arrange(.data$datetime)
  sorted_decimal_time <- sort(decimal_time)

  time_seq <- order(decimal_time)

  sorted_profile <-profile[time_seq]


  # Compute slopes as differences between consecutive melatonin concentration values
 # diff(profile)/diff(decimal_time)
  diff(sorted_profile)/diff(sorted_decimal_time)

}
