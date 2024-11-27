# Function to convert hms time to decimal hours
time_to_decimal <- function(hms_time, days_elapsed = 0) {
  # Extract hours, minutes, and seconds using lubridate
  hours <- hour(hms_time)   # Extract hours
  minutes <- minute(hms_time) # Extract minutes

  # Add days_elapsed * 24 to the hours and convert minutes to decimal
  total_hours <- (days_elapsed * 24) + hours + (minutes / 60)

  return(total_hours)
}

# Function to convert decimal hours back to hms format
# Function to convert decimal hours back to hms format, considering elapsed days
decimal_to_time <- function(decimal_time) {
  # Calculate the total elapsed days and hours from decimal time
  days <- floor(decimal_time / 24)            # Number of full days
  remaining_hours <- decimal_time %% 24       # Remainder hours within a single day

  # Extract hours and minutes
  hours <- floor(remaining_hours)              # Integer part for hours
  minutes <- (remaining_hours - hours) * 60    # Convert fractional part to minutes

  # Convert hours and minutes to seconds
  total_seconds <- (hours * 3600) + (minutes * 60)

  # Instead of using hms(), create a formatted string and use hms::as_hms() to parse it
  time_str <- sprintf("%02d:%02d:%02d", hours, floor(minutes), 0)

  # Create the hms object by parsing the formatted string
  hms_time <- as_hms(time_str)

  # Return elapsed days and the hms time
  return(list(days = days, time = hms_time))
}

# Convert POSIXct to Decimal Hours
posixct_to_decimal <- function(posix_times, origin_date = as.Date("2024-04-16")) {
  # Ensure input is POSIXct
  posix_times <- as.POSIXct(posix_times, tz = "UTC")

  # Calculate the number of elapsed days since origin
  days_elapsed <- as.numeric(as.Date(posix_times) - origin_date)
print(days_elapsed)
  # Extract the hours and minutes
  hours <- as.numeric(format(posix_times, "%H"))
  minutes <- as.numeric(format(posix_times, "%M"))

  # Calculate decimal hours for the current day
  decimal_time_today <- hours + (minutes / 60)

  # Total decimal hours, including elapsed days
  decimal_hours <- (days_elapsed * 24) + decimal_time_today
print(decimal_hours)
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

