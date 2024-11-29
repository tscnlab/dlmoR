# library(dplyr)
# library(ggplot2)
#
# # Function to convert POSIXct datetime to numeric (time difference in seconds from min(datetime))
# datetime_to_numeric <- function(datetime_col) {
#   # Find the minimum datetime (the reference point)
#   min_datetime <- min(datetime_col)
#
#   # Convert datetime to numeric (time difference in seconds)
#   return(as.numeric(difftime(datetime_col, min_datetime, units = "secs")))
# }
#
# # Function to convert numeric back to POSIXct datetime
# numeric_to_datetime <- function(numeric_value, datetime_reference) {
#   # Convert numeric value (seconds) back to datetime using the reference
#   return(as.POSIXct(numeric_value, origin = min(datetime_reference), tz = "UTC"))
# }
#
# # Get the corners of the parallelogram
# get_corners <- function(x0, y0, x1, y1, slope) {
#   if (slope == 0) {
#     lower_left <- c(x0, y0)
#     lower_right <- c(x1, y0)
#     upper_left <- c(x0, y1)
#     upper_right <- c(x1, y1)
#   } else {
#     height <- y1 - y0
#     delta_x <- height / slope
#     if (abs(slope) > 1e3) delta_x <- 0
#
#     lower_left <- c(x0, y0)
#     lower_right <- c(x1, y0)
#     upper_left <- c(x0 + delta_x, y1)
#     upper_right <- c(x1 + delta_x, y1)
#   }
#   return(list(lower_left, lower_right, upper_right, upper_left))
# }
#
# # Define the constraint functions for optimization
# constraints <- function(x, y, x0, x1, y0, y1, slope) {
#   corners <- get_corners(x0, y0, x1, y1, slope)
#   lower_left <- corners[[1]]
#   lower_right <- corners[[2]]
#   upper_right <- corners[[3]]
#   upper_left <- corners[[4]]
#
#   constraint_vals <- numeric()
#
#   for (i in seq_along(x)) {
#     xi <- x[i]
#     yi <- y[i]
#
#     y_lower <- if (xi <= lower_right[1]) y0 else y0 + slope * (xi - lower_right[1])
#     y_upper <- if (xi >= upper_left[1]) y1 else y1 + slope * (xi - upper_left[1])
#
#     x_lower <- min(lower_left[1], upper_left[1])
#     x_upper <- max(upper_right[1], lower_right[1])
#
#     constraint_vals <- c(constraint_vals, y_upper - yi)
#     constraint_vals <- c(constraint_vals, yi - y_lower)
#     constraint_vals <- c(constraint_vals, xi - x_lower)
#     constraint_vals <- c(constraint_vals, x_upper - xi)
#   }
#   return(constraint_vals)
# }
#
# # Objective function for optimization (parallelogram area)
# objective <- function(params, x, y, y0, y1) {
#   x0 <- params[1]
#   x1 <- params[2]
#   slope <- params[3]
#
#   corners <- get_corners(x0, y0, x1, y1, slope)
#
#   # Area calculation
#   v1 <- corners[[2]] - corners[[1]]
#   v2 <- corners[[4]] - corners[[1]]
#   area <- abs(v1[1] * v2[2] - v1[2] * v2[1])
#
#   # Apply penalty if constraints are violated
#   c_penalty <- min(constraints(x, y, x0, x1, y0, y1, slope))
#   if (c_penalty > 0) c_penalty <- 0
#
#   return(area + 1e3 * c_penalty^2)  # Added large penalty term for constraint violations
# }
#
# # Optimization function to fit the parallelogram
# optimize_parallelogram <- function(x, y) {
#   y0 <- min(y)
#   y1 <- max(y)
#
#   # Initial guess for optimization
#   x0_initial <- min(x)
#   x1_initial <- max(x)
#   slope_initial <- (y1 - y0) / (x1_initial - x0_initial)
#
#   initial_guess <- c(x0_initial, x1_initial, slope_initial)
#
#   # Optimize using L-BFGS-B method
#   result <- optim(
#     par = initial_guess,
#     fn = objective,
#     method = "L-BFGS-B",
#     lower = c(-Inf, -Inf, 0),
#     upper = c(Inf, Inf, Inf),
#     x = x,
#     y = y,
#     y0 = y0,
#     y1 = y1
#   )
#
#   # Return the optimized parameters
#   return(result$par)
# }
#
# # Function to apply the parallelogram truncation and convert back to datetime
# parallelogram_truncation <- function(profile_data) {
#   # Filter the profile_data to use only the rows where ascending == 1
#   profile_data_ascending <- profile_data %>%filter(ascending == 1)
#
#   # Convert datetime to numeric for optimization
#   x_values <- datetime_to_numeric(profile_data_ascending$datetime)
#   y_values <- profile_data_ascending$melatonin  # Use melatonin data (assuming it is in the 'melatonin' column)
#
#   # Get the optimized parameters for the parallelogram
#   params <- optimize_parallelogram(x_values, y_values)
#
#   # Extract the optimized x0 and x1 values (in numeric form)
#   optimized_x0 <- params[1]
#   optimized_x1 <- params[2]
#   slope <- params[3]
#
#
#   # Convert the optimized numeric x0 and x1 back to datetime using the reference (min of original datetime column)
#   trunc_datetime_0 <- numeric_to_datetime(optimized_x0, profile_data$datetime)
#   trunc_datetime_1 <- numeric_to_datetime(optimized_x1, profile_data$datetime)
#
#   # Print the result
#   print(trunc_datetime_0)
#   print(trunc_datetime_1)
#
#   # Return the result as a list with datetime values
#   return(list(trunc_datetime_0 = trunc_datetime_0, trunc_datetime_1 = trunc_datetime_1, slope))
# }
#
# # Example usage:
# # Assuming 'profile_data' is your tibble that contains 'datetime', 'melatonin', and 'ascending' columns:
# # Call the function to truncate and get the optimized time points
# # trunc_result <- parallelogram_truncation(profile_data)
# #
# # # Print the result
# # print(trunc_result)
# #
# # # Plot the original profile data with the optimized parallelogram on top
# # ggplot(profile_data, aes(x = datetime, y = melatonin)) +
# #   geom_line() +  # Plot the full melatonin profile
# #   geom_vline(xintercept = as.numeric(trunc_result$trunc_datetime_0), color = "red", linetype = "dashed") +
# #   geom_vline(xintercept = as.numeric(trunc_result$trunc_datetime_1), color = "red", linetype = "dashed") +
# #   labs(title = "Optimized Parallelogram Fit on Melatonin Profile",
# #        x = "Datetime", y = "Melatonin") +
# #   theme_minimal()
# Function to convert POSIXct datetime to numeric (time difference in seconds from min(datetime))
# Function to convert POSIXct datetime to numeric (time difference in hours from min(datetime))
datetime_to_numeric <- function(datetime_col, datetime_reference) {
  # Find the minimum datetime (the reference point)
  min_datetime <- min(datetime_reference)

  # Convert datetime to numeric (time difference in hours)
  return(as.numeric(difftime(datetime_col, min_datetime, units = "hours")))
}

# Function to convert numeric back to POSIXct datetime
numeric_to_datetime <- function(numeric_value, datetime_reference) {
  # Convert numeric value (hours) back to datetime using the reference
  min_datetime <- min(datetime_reference)
  return(min_datetime + lubridate::dhours(numeric_value))
}

# omg when was i last using this? 27.11.2024

# # Get the corners of the parallelogram
# get_corners <- function(x0, y0, x1, y1, slope) {
#   if (slope == 0) {
#     lower_left <- c(x0, y0)
#     lower_right <- c(x1, y0)
#     upper_left <- c(x0, y1)
#     upper_right <- c(x1, y1)
#   } else {
#     height <- y1 - y0
#     delta_x <- height / slope
#     if (abs(slope) > 1e3) delta_x <- 0
#
#     lower_left <- c(x0, y0)
#     lower_right <- c(x1, y0)
#     upper_left <- c(x0 + delta_x, y1)
#     upper_right <- c(x1 + delta_x, y1)
#   }
#   return(list(lower_left, lower_right, upper_right, upper_left))
# }
#
# # Define the constraint functions for optimization
# constraints <- function(x, y, x0, x1, y0, y1, slope) {
#   corners <- get_corners(x0, y0, x1, y1, slope)
#   lower_left <- corners[[1]]
#   lower_right <- corners[[2]]
#   upper_right <- corners[[3]]
#   upper_left <- corners[[4]]
#
#   constraint_vals <- numeric()
#
#   for (i in seq_along(x)) {
#     xi <- x[i]
#     yi <- y[i]
#
#     y_lower <- if (xi <= lower_right[1]) y0 else y0 + slope * (xi - lower_right[1])
#     y_upper <- if (xi >= upper_left[1]) y1 else y1 + slope * (xi - upper_left[1])
#
#     x_lower <- min(lower_left[1], upper_left[1])
#     x_upper <- max(upper_right[1], lower_right[1])
#
#     constraint_vals <- c(constraint_vals, y_upper - yi)
#     constraint_vals <- c(constraint_vals, yi - y_lower)
#     constraint_vals <- c(constraint_vals, xi - x_lower)
#     constraint_vals <- c(constraint_vals, x_upper - xi)
#   }
#   return(constraint_vals)
# }
#
# # Objective function for optimization (parallelogram area)
# objective <- function(params, x, y, y0, y1) {
#   x0 <- params[1]
#   x1 <- params[2]
#   slope <- params[3]
#
#   corners <- get_corners(x0, y0, x1, y1, slope)
#
#   # Area calculation
#   v1 <- corners[[2]] - corners[[1]]
#   v2 <- corners[[4]] - corners[[1]]
#   area <- abs(v1[1] * v2[2] - v1[2] * v2[1])
#
#   # Apply penalty if constraints are violated
#   c_penalty <- min(constraints(x, y, x0, x1, y0, y1, slope))
#   if (c_penalty > 0) c_penalty <- 0
#
#   return(area + 1e3 * c_penalty^2)  # Added large penalty term for constraint violations
# }
#
# # Optimization function to fit the parallelogram
# optimize_parallelogram <- function(x, y) {
#   y0 <- min(y)
#   y1 <- max(y)
#
#   # Initial guess for optimization
#   x0_initial <- min(x)
#   x1_initial <- max(x)
#   slope_initial <- (y1 - y0) / (x1_initial - x0_initial)
#
#   initial_guess <- c(x0_initial, x1_initial, slope_initial)
#
#   # Optimize using L-BFGS-B method
#   result <- optim(
#     par = initial_guess,
#     fn = objective,
#     method = "L-BFGS-B",
#     lower = c(-Inf, -Inf, 0),
#     upper = c(Inf, Inf, Inf),
#     x = x,
#     y = y,
#     y0 = y0,
#     y1 = y1
#   )
#
#   # Return the optimized parameters
#   return(result$par)
# }

# Function to apply the parallelogram truncation and convert back to datetime
parallelogram_truncation <- function(profile_data) {
  # Filter the profile_data to use only the rows where ascending == 1
  profile_data_ascending <- profile_data %>% filter(ascending == 1)
  print(profile_data_ascending)
  # Convert datetime to numeric for optimization
  # x_values <- datetime_to_numeric(profile_data_ascending$datetime, profile_data$datetime)
  # x_values <- datetime_to_numeric(profile_data_ascending$datetime, profile_data$datetime)
  x_values <- posixct_to_decimal(profile_data_ascending$datetime)

  print("x_values")
  print(x_values)
  y_values <- profile_data_ascending$melatonin  # Use melatonin data (assuming it is in the 'melatonin' column)
  print("y_values")
  print(y_values)

  # Get the optimized parameters for the parallelogram
  params <- optimize_parallelogram(x_values, y_values)
  print("optimized parameters")
  print(params)
  # Extract the optimized x0 and x1 values (in numeric form)
  optimized_x0 <- params[1]
  optimized_x1 <- params[2]
  slope <- params[3]

  # Convert the optimized numeric x0 and x1 back to datetime using the reference (min of original datetime column)
  # trunc_datetime_0 <- numeric_to_datetime(optimized_x0, profile_data$datetime)
  # trunc_datetime_1 <- numeric_to_datetime(optimized_x1, profile_data$datetime)
    trunc_datetime_0 <- decimal_to_posixct(optimized_x0)
    trunc_datetime_1 <- decimal_to_posixct(optimized_x1)
  # Print the result
  print(trunc_datetime_0)
  print(trunc_datetime_1)

  # Return the result as a list with datetime values
  return(list(trunc_datetime_0 = trunc_datetime_0, trunc_datetime_1 = trunc_datetime_1, pll_slope = slope))
}
