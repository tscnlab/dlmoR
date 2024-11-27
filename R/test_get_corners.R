# Minimal Test: Check get_corners with datetime
test_get_corners <- function() {
  # Create some test data for x0, y0, x1, y1, slope
  x0 <- as.POSIXct("2024-04-16 20:23:20", tz = "UTC")  # POSIXct datetime
  x1 <- as.POSIXct("2024-04-16 23:10:00", tz = "UTC")  # POSIXct datetime
  y0 <- 10  # Example melatonin value
  y1 <- 100  # Example melatonin value
  slope <- 0.005  # Example slope

  # Call get_corners
  corners <- get_corners(x0, y0, x1, y1, slope)

  # Check the results
  print(corners)  # Check the output

  # Convert the corners back to POSIXct for plotting
  lower_left_time <- as.POSIXct(corners[[1]][1], origin = "1970-01-01", tz = "UTC")
  lower_right_time <- as.POSIXct(corners[[2]][1], origin = "1970-01-01", tz = "UTC")
  upper_left_time <- as.POSIXct(corners[[4]][1], origin = "1970-01-01", tz = "UTC")
  upper_right_time <- as.POSIXct(corners[[3]][1], origin = "1970-01-01", tz = "UTC")

  print(lower_left_time)
  print(lower_right_time)
  print(upper_left_time)
  print(upper_right_time)

  return(corners)
}

# Run the test function
#test_get_corners()
