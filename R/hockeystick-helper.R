#' Calculate Slopes Between Data Points
#'
#' This function calculates the slopes (rate of change) between consecutive data points
#' in a melatonin profile based on melatonin concentrations and their corresponding
#' time points.
#'
#' @param profile Numeric vector representing melatonin concentrations at each time point.
#' @param decimal_time Numeric vector representing the time points in decimal format (e.g., fractional hours or days).
#' @return A numeric vector of slopes calculated as the change in melatonin concentration
#' divided by the change in time between consecutive points.
#'
#' @details
#' The function sorts the input data by time to ensure proper calculation of slopes.
#' The slopes are computed as:
#' \deqn{slope = (profile[i+1] - profile[i]) / (decimal_time[i+1] - decimal_time[i])}.
#'
#' @examples
#' profile <- c(1.2, 1.4, 1.5, 1.7, 2.0)
#' decimal_time <- c(0, 0.25, 0.5, 0.75, 1)  # Time in fractional hours
#' slopes <- calculate_slopes(profile, decimal_time)
#' print(slopes)
#' @export
calculate_slopes <- function(profile, decimal_time) {
  # Ensure the time points are sorted
  sorted_decimal_time <- sort(decimal_time)

  # Reorder the profile data to match the sorted time points
  time_seq <- order(decimal_time)  # Get the order of indices for sorting
  sorted_profile <- profile[time_seq]

  # Compute slopes as the rate of change in melatonin concentrations over time
  slopes <- diff(sorted_profile) / diff(sorted_decimal_time)

  # Return the calculated slopes
  return(slopes)
}

#' Convert POSIXct Timestamps to Decimal Hours
#'
#' This function converts a vector of POSIXct timestamps into decimal hours,
#' relative to the first timestamp in `profile_datetime`.
#'
#' @param posix_times A vector of POSIXct timestamps.
#' @param profile_datetime A vector containing POSIXct timestamps, where the first element serves as the reference origin.
#' @return A numeric vector representing decimal hours since the origin.
#' @export
posixct_to_decimal <- function(posix_times, profile_datetime) {
  # Ensure input is in POSIXct format
  posix_times <- as.POSIXct(posix_times, tz = "UTC")

  # Extract the first timestamp as the reference origin
  posix_origin <- profile_datetime[1]
  origin_date <- as.Date(posix_origin)  # Convert to Date format

  # Compute elapsed days from the origin
  days_elapsed <- as.numeric(as.Date(posix_times) - origin_date)

  # Extract hours, minutes, and seconds from the timestamps
  hours <- as.numeric(format(posix_times, "%H"))
  minutes <- as.numeric(format(posix_times, "%M"))
  seconds <- as.numeric(format(posix_times, "%S"))

  # Convert time of day to decimal hours
  decimal_time_today <- hours + (minutes / 60) + (seconds / 3600)

  # Compute total decimal hours, including days elapsed
  decimal_hours <- (days_elapsed * 24) + decimal_time_today

  return(decimal_hours)
}


#' Convert Decimal Hours to POSIXct Timestamps
#'
#' This function converts a numeric vector of decimal hours back into POSIXct timestamps,
#' using the first timestamp in `profile_datetime` as the origin.
#'
#' @param decimal_hours A numeric vector representing time in decimal hours.
#' @param profile_datetime A vector of POSIXct timestamps where the first element serves as the reference origin.
#' @param tz A string specifying the time zone (default is "UTC").
#' @return A vector of POSIXct timestamps.
#' @export
decimal_to_posixct <- function(decimal_hours, profile_datetime, tz = "UTC") {
  # Compute the number of full days elapsed
  days_elapsed <- floor(decimal_hours / 24)

  # Compute the remaining hours within the last day
  remaining_hours <- decimal_hours %% 24

  # Extract hours, minutes, and seconds from the decimal hours
  hour <- floor(remaining_hours)
  minute_fraction <- (remaining_hours - hour) * 60
  minute <- floor(minute_fraction)
  second <- round((minute_fraction - minute) * 60)  # Round to nearest second

  # Extract the origin date from the profile_datetime
  posix_origin <- profile_datetime[1]
  origin_date <- as.Date(posix_origin)  # Ensure it's just the date

  # Compute full POSIXct timestamp
  posix_time <- as.POSIXct(origin_date, tz = tz) +
    days_elapsed * 86400 +  # seconds in a day
    hour * 3600 +
    minute * 60 +
    second

  return(posix_time)
}

