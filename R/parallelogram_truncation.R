#' Apply Parallelogram Truncation and Convert Back to Datetime
#'
#' This function fits a parallelogram to the ascending segment of melatonin profile data,
#' calculates its corners, and determines whether it meets predefined constraints.
#' If the parallelogram does not satisfy the slope ratio criteria, a warning is issued.
#'
#' @param profile_data A tibble containing melatonin concentration data with datetime values.
#' @return A list containing:
#'   - `pll_datetime_0`: The leftmost x-coordinate (datetime) of the parallelogram.
#'   - `pll_datetime_1`: The rightmost x-coordinate (datetime) of the parallelogram.
#'   - `pll_slope`: The optimized slope of the parallelogram.
#'   - `corners`: The four corner points of the parallelogram.
#'   - `flag`: Logical indicating whether the parallelogram violates constraints.
#'
#' @seealso \code{\link{optimize_parallelogram}}, \code{\link{get_corners}}
#' @export
parallelogram_fit <- function(profile_data) {
  flag <- FALSE  # Initialize flag to track constraint violations

  # Filter the profile_data to extract only ascending data points
  profile_data_ascending <- profile_data %>% dplyr::filter(ascending == 1)

  # Convert datetime to numeric (decimal time) for optimization
  x_values <- posixct_to_decimal(profile_data_ascending$datetime, profile_data$datetime[3])
  y_values <- profile_data_ascending$melatonin

  # Identify the first ascending index in the profile data
  first_ascending_index <- which(profile_data$ascending == 1)[1]
  # Identify the last point before the ascending segment
  last_point_before_ascending_index <- first_ascending_index - 1
  lpba <- profile_data[last_point_before_ascending_index, ]  # Extract last point before ascending
  lpba_dec <- posixct_to_decimal(lpba$datetime, profile_data$datetime[3])

  # Append the last point before ascending to the data used for optimization
  x_pll <- c(lpba_dec, x_values)
  y_pll <- c(lpba$melatonin, y_values)

  # Determine vertical boundaries (minimum and maximum y-values)
  y0 <- min(y_pll)
  y1 <- max(y_pll)

  # Optimize the parallelogram based on x and y values
  params <- optimize_parallelogram(x_pll, y_pll)

  # Extract optimized x0, x1, and slope values
  x0_numeric <- params[1]
  x1_numeric <- params[2]
  slope <- params[3]

  # Compute corners of the parallelogram
  corners <- get_corners(x0_numeric, y0, x1_numeric, y1, slope)

  # Convert the list of corner points into a tibble for further analysis
  corners_tibble <- purrr::map_dfr(
    corners,
    ~ tibble::tibble(x = .[1], y = .[2]),
    .id = "corner"
  )

  # Compute the ratio of lateral side slope to diagonal slope
  diagonal_slopes <- corners_tibble %>%
    dplyr::summarize(
      slope_diag1 = (y[which(corner == "ur")] - y[which(corner == "lr")]) /  # Right lateral side
        (x[which(corner == "ur")] - x[which(corner == "lr")]),
      slope_diag2 = (y[which(corner == "ur")] - y[which(corner == "ll")]) /  # Left diagonal
        (x[which(corner == "ur")] - x[which(corner == "ll")])
    ) %>%
    dplyr::mutate(ratio = (slope_diag1 / slope_diag2))

  # Check for constraint violations
  if (diagonal_slopes$ratio < 0) {
    warning("The ratio of the slopes is negative.")
    flag <- TRUE
  }
  if (diagonal_slopes$ratio < 0.48) {  # Threshold constraint for parallelogram shape
    warning("The diagonal slope ratio is less than 0.48.")
    flag <- TRUE
  }

  # Convert numeric x0 and x1 back to datetime format
  pll_datetime_0 <- decimal_to_posixct(x0_numeric, profile_data$datetime[3])
  pll_datetime_1 <- decimal_to_posixct(x1_numeric, profile_data$datetime[3])

  # Return results as a list
  return(list(
    pll_datetime_0 = pll_datetime_0,
    pll_datetime_1 = pll_datetime_1,
    pll_slope = slope,
    corners = corners,
    flag = flag
  ))
}
