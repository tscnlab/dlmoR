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
    lower_left <- c(x0, y0)
    lower_right <- c(x1, y0)
    upper_left <- c(x0 + delta_x, y1)
    upper_right <- c(x1 + delta_x, y1)
  }
  return(list(ll = lower_left, lr =lower_right, ur = upper_right, ul = upper_left))
}


# Helper Function for the Constraints
constraints <- function(params, x, y, y0, y1) {
  x0 <- params[1]
  x1 <- params[2]
  slope <- params[3]

  corners <- get_corners(x0, y0, x1, y1, slope)
  lower_left <- corners[[1]]
  lower_right <- corners[[2]]
  upper_right <- corners[[3]]
  upper_left <- corners[[4]]

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

  x_lower <- min(lower_left[1], upper_left[1])
  x_upper <- max(lower_right[1], upper_right[1])

  # Compute constraint violations
  c(
    y_upper - y,       # Point is below the upper edge
    y - y_lower,       # Point is above the lower edge
    x - x_lower,       # Point is right of the left edge
    x_upper - x        # Point is left of the right edge
  )

}

# Helper Function for the Objective
objective <- function(params, x, y, y0, y1) {
  x0 <- params[1]
  x1 <- params[2]
  slope <- params[3]

  corners <- get_corners(x0, y0, x1, y1, slope)
  width <- abs(corners$lr[1] - corners$ll[1])
  height <- abs(corners$ul[2] - corners$ll[2])
  area <- width * height
  c_penalty <- min(constraints(params, x, y, y0, y1))
  if (c_penalty >= 0) {
    c_penalty <- 0
  }
  if(corners$lr[1] - corners$ll[1] > 0){
    c_left_right <-0 # left is smaller than right (i.e., left is left of right)
  } else{
    c_left_right <-  (corners$lr[1] - corners$ll[1]) ** 2
  }
  return(area + 1e3 * log(1 + c_penalty^2) + 1e3 * log(c_left_right+1))
}

# Main Function to Optimize Parallelogram
optimize_parallelogram <- function(x, y) {
  y0 <- min(y)
  y1 <- max(y)

  x0_initial <- x[1]
  x_slope_2 <- x[2]
  x1_initial <- tail(x, n=1)

  slope_initial <- (y1 - y0) / (x_slope_2 - x0_initial)
  initial_guess <- c(x0_initial*0.8, x1_initial*1.2, slope_initial)

  result <- stats::optim(
    par = initial_guess,
    fn = function(params) objective(params, x, y, y0, y1),
    method = "SANN",
    control=list(maxit=100000)
  )


    x0 <- result$par[1]
  x1 <- result$par[2]
  slope <- result$par[3]

  corners <- get_corners(x0, y0, x1, y1, slope)
  #print(corners)

  result$par
}


# optimize_parallelogram(
#   c(18.53333, 19.28333, 20.78333),
#   c(0.100, 3.174, 20.109)
# )


# optimize_parallelogram(
#   c(17, 18.5, 19.5),
#   c(0.100, 0.1, 12)
# )

# optimize_parallelogram(
#   c(17, 18, 19.5, 20, 21, 22, 23),
#   c(1.1, 0.1, 7.5, 13, 14, 16, 19)
# )

# optimize_parallelogram(
#   c(18.56, 20.05),
#   c(1.84, 9.203)
# )

