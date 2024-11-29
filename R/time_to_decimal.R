# Convert POSIXct to Decimal Hours
# posixct_to_decimal <- function(posix_times, origin_date = as.Date("2024-04-16")) {
posixct_to_decimal <- function(posix_times, profile_datetime) {

  # Ensure input is POSIXct
  posix_times <- as.POSIXct(posix_times, tz = "UTC")

  # Extract reference origin date
  posix_origin <-profile_datatime[1]
  origin_date<-as.Date(posix_origin)

  # Calculate the number of elapsed days since origin
  days_elapsed <- as.numeric(as.Date(posix_times) - origin_date)

  # Extract the hours and minutes
  hours <- as.numeric(format(posix_times, "%H"))
  minutes <- as.numeric(format(posix_times, "%M"))

  # Calculate decimal hours for the current day
  decimal_time_today <- hours + (minutes / 60)

  # Total decimal hours, including elapsed days
  decimal_hours <- (days_elapsed * 24) + decimal_time_today

  return(decimal_time_today)
}


# Convert Decimal Hours to POSIXct
decimal_to_posixct <- function(decimal_hours, origin_date = as.Date("2024-04-16"), tz = "UTC") {
  # Calculate the number of full days elapsed
  days_elapsed <- floor(decimal_hours / 24)

  # Remaining decimal hours within the current day
  remaining_hours <- decimal_hours %% 24

  # Extract hours and minutes
  hour <- floor(remaining_hours)
  minute <- (remaining_hours - hour) * 60

  # Create POSIXct for the origin date + elapsed days + time of day
  posix_time <- as.POSIXct(origin_date, tz = tz) + days_elapsed * 86400 + hour * 3600 + minute * 60

  return(posix_time)
}

