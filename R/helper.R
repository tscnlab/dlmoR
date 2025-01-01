# calculate slopes between data points
calculate_slopes <- function(profile, decimal_time) {
  # Ensure datetime is sorted
  # Ensure datetime is sorted
  # profile_data <- profile_data %>% dplyr::arrange(.data$datetime)
  print(decimal_time)
  sorted_decimal_time <- sort(decimal_time)
  print(sorted_decimal_time)
  time_seq <- order(decimal_time)
  print(profile)
  sorted_profile <-profile[time_seq]
  print(sorted_profile)

  # Compute slopes as differences between consecutive melatonin concentration values
 # diff(profile)/diff(decimal_time)
  diff(sorted_profile)/diff(sorted_decimal_time)
  print(diff(sorted_profile)/diff(sorted_decimal_time))
}
