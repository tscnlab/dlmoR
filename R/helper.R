# calculate slopes between data points
calculate_slopes <- function(profile) {
  # Compute slopes as differences between consecutive melatonin concentration values
  diff(profile)
}
