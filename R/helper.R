# calculate slopes between data points
calculate_slopes <- function(profile, decimal_time) {
  # Compute slopes as differences between consecutive melatonin concentration values
  diff(profile)/diff(decimal_time)
}
