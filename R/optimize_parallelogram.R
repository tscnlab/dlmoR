# Helper Function to Get Corners of the Parallelogram
get_corners <- function(x0, y0, x1, y1, slope) {
  if (slope == 0) {
    lower_left <- c(x0, y0)
    lower_right <- c(x1, y0)
    upper_left <- c(x0, y1)
    upper_right <- c(x1, y1)
  } else {
    height <- y1 - y0
    delta_x <- height / slope

    if (abs(slope) > 1e3) {
      delta_x <- 0
    }

    lower_left <- c(x0, y0)
    lower_right <- c(x1, y0)
    upper_left <- c(x0 + delta_x, y1)
    upper_right <- c(x1 + delta_x, y1)
  }
  return(list(lower_left, lower_right, upper_right, upper_left))
}


# Helper Function for the Constraints
constraints <- function(params, x, y, y0, y1) {
  x0 <- params[1]
  x1 <- params[2]
  slope <- params[3]
  constraint_vals <- numeric()

  corners <- get_corners(x0, y0, x1, y1, slope)
  lower_left <- corners[[1]]
  lower_right <- corners[[2]]
  upper_right <- corners[[3]]
  upper_left <- corners[[4]]

  # Loop over x and y values to compute all constraints
  for (i in seq_along(x)) {
    xi <- x[i]
    yi <- y[i]

    y_lower <- if(xi <= lower_right[1]) y0 else y0 + slope * (xi - lower_right[1])
    y_upper <- if(xi >= upper_left[1]) y1 else y1 + slope * (xi - upper_left[1])

    x_lower <- min(lower_left[1], upper_left[1])
    x_upper <- max(upper_right[1], lower_right[1])

    # Append all constraints
    constraint_vals <- c(
      constraint_vals,
      y_upper - yi,
      yi - y_lower,
      xi - x_lower,
      x_upper - xi
    )
  }
  constraint_vals
}

# Helper Function for the Objective
objective <- function(params, x, y, y0, y1) {
  x0 <- params[1]
  x1 <- params[2]
  slope <- params[3]

  corners <- get_corners(x0, y0, x1, y1, slope)
  v1 <- corners[[2]] - corners[[1]]
  v2 <- corners[[4]] - corners[[1]]
  area <- abs(v1[1] * v2[2] - v1[2] * v2[1])

  c_penalty <- min(constraints(params, x, y, y0, y1))
  if (c_penalty > 0) {
    c_penalty <- 0
  }

  return(area + 1e3 * c_penalty^2)
}

# Main Function to Optimize Parallelogram
optimize_parallelogram <- function(x, y) {
  y0 <- min(y)
  y1 <- max(y)
  x0_initial <- min(x)*0.9
  x1_initial <- max(x)*1.1
  slope_initial <- (y1 - y0) / (x1_initial - x0_initial)

  initial_guess <- c(x0_initial, x1_initial, slope_initial)

  result <- stats::optim(
    par = initial_guess,
    fn = function(params) objective(params, x, y, y0, y1),
    method = "L-BFGS-B",
    lower = c(-Inf, -Inf, 0),
    upper = c(Inf, Inf, Inf),
    control = list(maxit=1000, factr=1e8, trace=3)
  )

  result$par

}
