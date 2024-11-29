# Function to apply the parallelogram truncation and convert back to datetime
parallelogram_truncation <- function(profile_data) {

  # Filter the profile_data to use only the rows where ascending == 1
  profile_data_ascending <- profile_data %>% dplyr::filter(ascending == 1)

  # Convert datetime to numeric for optimization
  x_values <- posixct_to_decimal(profile_data_ascending$datetime, profile_data$datetime)
  y_values <- profile_data_ascending$melatonin

  # Get the optimized parameters for the parallelogram
  params <- optimize_parallelogram(x_values, y_values)

  # Extract the optimized x0 and x1 values (in numeric form)
  optimized_x0 <- params[1]
  optimized_x1 <- params[2]
  slope <- params[3]

  # Convert the optimized numeric x0 and x1 back to datetime
    pll_datetime_0 <- decimal_to_posixct(optimized_x0, profile_data$datetime)
    pll_datetime_1 <- decimal_to_posixct(optimized_x1, profile_data$datetime)

  # Return the result as a list with datetime values
  return(list(trunc_datetime_0 = pll_datetime_0, pll_datetime_1 = pll_datetime_1, pll_slope = slope))
}
