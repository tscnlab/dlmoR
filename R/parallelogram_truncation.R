# Function to apply the parallelogram truncation and convert back to datetime
parallelogram_fit <- function(profile_data) {
  flag<-FALSE

  # Filter the profile_data to use only the rows where ascending == 1
  profile_data_ascending <- profile_data %>% dplyr::filter(ascending == 1)


  # Convert datetime to numeric for optimization
  x_values <- posixct_to_decimal(profile_data_ascending$datetime, profile_data$datetime)
  y_values <- profile_data_ascending$melatonin

  # Find the row index of the first point in 'ascending' within 'profile_data'
  first_ascending_index <- which(profile_data$ascending == 1)[1]
  # Get the index of the last point before 'ascending'
  last_point_before_ascending_index <- first_ascending_index - 1
  # Filter the last point before the ascending subset
  lpba <- profile_data[last_point_before_ascending_index, ]
  lpba_dec <- posixct_to_decimal(lpba$datetime, profile_data$datetime)

  x_pll <- c(lpba_dec, x_values)
  y_pll <- c(lpba$melatonin, y_values)

  #print("xll")
  #print(x_pll)
  #print("yll")
  #print(y_pll)
 # y0 <- min(profile_data_ascending$melatonin)
#  y1 <- max(profile_data_ascending$melatonin)
  y0 <- min(y_pll)
  y1 <- max(y_pll)
  # Get the optimized parameters for the parallelogram
  # params <- optimize_parallelogram(x_values, y_values)
  params <- optimize_parallelogram(x_pll, y_pll)


  # Extract the optimized x0 and x1 values (in numeric form)
  x0_numeric <- params[1]
  x1_numeric <- params[2]
  slope <- params[3]

  # Get corners of the parallelogram
  corners <- get_corners(x0_numeric, y0, x1_numeric, y1, slope)

  corners_tibble <- purrr::map_dfr(
    corners,
    ~ tibble::tibble(x = .[1], y = .[2]),
    .id = "corner"
  )


  # corners_tibble <- corners_tibble %>%
  #   dplyr::mutate(x = dplyr::if_else(corner == "ur", 26, x))
  # Calculate slopes of the diagonals
  diagonal_slopes <- corners_tibble %>%
    dplyr::summarize(
      slope_diag1 = (y[which(corner == "ur")] - y[which(corner == "ll")]) /
        (x[which(corner == "ur")] - x[which(corner == "ll")]),
      slope_diag2 = (y[which(corner == "ul")] - y[which(corner == "lr")]) /
        (x[which(corner == "ul")] - x[which(corner == "lr")])
    ) %>%
    #dplyr::mutate(ratio = abs(slope_diag1 / slope_diag2))
    dplyr::mutate(ratio = (slope_diag1 / slope_diag2))
    #dplyr::mutate(ratio = (slope_diag2 / slope_diag1))


  # print("diag slope ratio")
  # #print("inverted diag slope ratio")
  #
  # print(diagonal_slopes$ratio)
  # print("corners")
  # print(corners_tibble)
  # Add warnings
  if (diagonal_slopes$ratio < 0) {
    warning("The ratio of the slopes is negative.")
    flag <- TRUE
  }
  if(diagonal_slopes$ratio < 0.5){
    warning("the diagonal slope ratio is less than 1/2")
    flag <- TRUE
  }

  #print(diagonal_slopes$ratio)
  # Convert the optimized numeric x0 and x1 back to datetime
    pll_datetime_0 <- decimal_to_posixct(x0_numeric, profile_data$datetime)
    pll_datetime_1 <- decimal_to_posixct(x1_numeric, profile_data$datetime)

  # Return the result as a list with datetime values
  return(list(pll_datetime_0 = pll_datetime_0, pll_datetime_1 = pll_datetime_1, pll_slope = slope, corners = corners, flag = flag))
}
