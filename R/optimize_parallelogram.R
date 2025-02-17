#' Compute the Four Corner Points of a Parallelogram
#'
#' Given the bottom-left and top-right coordinates, this function calculates
#' the four corners of a parallelogram with a specified slope.
#'
#' @param x0 Numeric. The x-coordinate of the bottom-left corner.
#' @param y0 Numeric. The y-coordinate of the bottom-left corner.
#' @param x1 Numeric. The x-coordinate of the bottom-right corner.
#' @param y1 Numeric. The y-coordinate of the top-right corner.
#' @param slope Numeric. The slope of the parallelogram's diagonal edges.
#' @return A list containing the coordinates of the four corners.
#' @export
get_corners <- function(x0, y0, x1, y1, slope) {
  if (slope == 0) {
    # If slope is zero, the parallelogram is a rectangle.
    lower_left <- c(x0, y0)
    lower_right <- c(x1, y0)
    upper_left <- c(x0, y1)
    upper_right <- c(x1, y1)
  } else {
    # Compute horizontal shift based on slope and height difference
    height <- y1 - y0
    delta_x <- height / slope
    lower_left <- c(x0, y0)
    lower_right <- c(x1, y0)
    upper_left <- c(x0 + delta_x, y1)
    upper_right <- c(x1 + delta_x, y1)
  }
  return(list(ll = lower_left, lr = lower_right, ur = upper_right, ul = upper_left))
}

#' Compute Constraint Violations for Parallelogram Fit
#'
#' Given a set of parameters, this function computes constraint violations,
#' ensuring that all points remain inside the parallelogram.
#'
#' @param params Numeric vector of optimization parameters: x0, x1, slope.
#' @param x Numeric vector of x-coordinates for data points.
#' @param y Numeric vector of y-coordinates for data points.
#' @param y0 Numeric. The minimum y-value in the dataset.
#' @param y1 Numeric. The maximum y-value in the dataset.
#' @return A numeric vector of constraint violations.
#' @export
constraints <- function(params, x, y, y0, y1) {
  x0 <- params[1]
  x1 <- params[2]
  slope <- params[3]

  # Get parallelogram corners based on computed parameters
  corners <- get_corners(x0, y0, x1, y1, slope)
  lower_left <- corners$ll
  lower_right <- corners$lr
  upper_right <- corners$ur
  upper_left <- corners$ul

  # Compute vertical boundary constraints
  y_lower <- ifelse(
    x <= lower_right[1],
    y0,
    y0 + slope * (x - lower_right[1])
  )
  y_upper <- ifelse(
    x >= upper_left[1],
    y1,
    y1 + slope * (x - upper_left[1])
  )

  # Compute horizontal boundary constraints
  x_lower <- min(lower_left[1], upper_left[1])
  x_upper <- max(lower_right[1], upper_right[1])

  # Return constraint violations
  c(
    y_upper - y,       # Points should not be above upper edge
    y - y_lower,       # Points should not be below lower edge
    x - x_lower,       # Points should not be to the left of left edge
    x_upper - x        # Points should not be to the right of right edge
  )
}

#' Objective Function for Parallelogram Optimization
#'
#' This function defines the objective to minimize when fitting a parallelogram.
#' The objective is to find the smallest parallelogram that includes all points.
#'
#' @param params Numeric vector of optimization parameters: x0, x1, slope.
#' @param x Numeric vector of x-coordinates for data points.
#' @param y Numeric vector of y-coordinates for data points.
#' @param y0 Numeric. The minimum y-value in the dataset.
#' @param y1 Numeric. The maximum y-value in the dataset.
#' @return The computed objective value (area plus penalties for constraint violations).
#' @export
objective <- function(params, x, y, y0, y1) {
  x0 <- params[1]
  x1 <- params[2]
  slope <- params[3]

  # Compute parallelogram corners
  corners <- get_corners(x0, y0, x1, y1, slope)
  width <- abs(corners$lr[1] - corners$ll[1])  # Horizontal width
  height <- abs(corners$ul[2] - corners$ll[2]) # Vertical height
  area <- width * height

  # Compute constraint penalties
  c_penalty <- min(constraints(params, x, y, y0, y1))
  if (c_penalty >= 0) {
    c_penalty <- 0
  }
  if (corners$lr[1] - corners$ll[1] > 0) {
    c_left_right <- 0 # Ensure left is smaller than right
  } else {
    c_left_right <- (corners$lr[1] - corners$ll[1])^2
  }

  # Return total cost (area + penalties)
  return(area + 1e3 * log(1 + c_penalty^2) + 1e3 * log(c_left_right + 1))
}


#' Main Function to Optimize Parallelogram
#'
#' Optimize Parallelogram to Fit Input Data
#'
#' This function fits the smallest possible parallelogram that fully encloses all input points.
#' It minimizes the area while ensuring that all points remain within the defined region.
#'
#' @param x A numeric vector representing x-coordinates (e.g., time in decimal hours).
#' @param y A numeric vector representing y-coordinates (e.g., melatonin concentration).
#' @return A list containing optimized parameters for the parallelogram, including corner points.
#' @examples
#' optimize_parallelogram(c(18.56, 20.05), c(1.84, 9.203))
#' @export

optimize_parallelogram <- function(x, y) {
  # Identify the minimum and maximum y-values
  y0 <- min(y)
  y1 <- max(y)

  # Initialize parameters for optimization
  x0_initial <- x[1]
  x1_initial <- tail(x, n=1)
  slope_initial <- (y1 - y0) / (x1_initial - x0_initial)
  delta_y <- y1 - y0
  x_right_initial <- x1_initial - delta_y / slope_initial
  initial_guess <- c(x0_initial * 0.8, x_right_initial * 1.1, slope_initial)

  # Perform optimization to minimize area while including all points
  result <- stats::optim(
    par = initial_guess,
    fn = function(params) objective(params, x, y, y0, y1),
    method = "SANN",
    control = list(maxit = 200000) # Simulated Annealing for global optimization
  )

  # Extract optimized parameters
  x0 <- result$par[1]
  x1 <- result$par[2]
  slope <- result$par[3]

  return(result$par)
}
