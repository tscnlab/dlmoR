# Load required libraries
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

# Function to convert POSIXct times to decimal hours
posixct_to_decimal <- function(posix_times, profile_datetime) {
  posix_times <- as.POSIXct(posix_times, tz = "UTC")
  posix_origin <- profile_datetime[1]
  origin_date <- as.Date(posix_origin)

  days_elapsed <- as.numeric(as.Date(posix_times) - origin_date)
  hours <- as.numeric(format(posix_times, "%H"))
  minutes <- as.numeric(format(posix_times, "%M"))

  decimal_time_today <- hours + (minutes / 60)
  decimal_hours <- (days_elapsed * 24) + decimal_time_today
  return(decimal_hours)
}

# Define a function to create a grid of points within a region of interest
make_grid <- function(roi, step_x, step_y) {
  xmin <- roi$x[1]
  xmax <- roi$x[2]
  ymin <- roi$y[1]
  ymax <- roi$y[2]

  x_seq <- seq(from = xmin, to = xmax, by = step_x)
  y_seq <- seq(from = ymin, to = ymax, by = step_y)

  grid_points <- expand.grid(x = x_seq, y = y_seq)
  return(grid_points)
}

# Define a nonlinear constraint function for parabolic fitting
nl_constraints <- function(params, poi_x, poi_y, x) {
  a <- params[1]
  b <- params[2]
  c <- params[3]

  # Constraint 1: Parabola must pass through the POI
  poi_diff <- a * poi_x^2 + b * poi_x + c - poi_y

  # Constraint 2: Ensure positive slope (dy/dx > 0) for all points to the right of the POI
  slope_values <- 2 * a * x + b
  slope_min <- min(slope_values)

  if (slope_min >= 0) {
    slope_min <- 0
  }

  return(list(pass_poi = poi_diff^2, pos_grad = slope_min^2))
}

# Define the objective function for line or parabola fitting
objective_function <- function(params, x, y, poi, fit_type) {
  poi_x <- poi$x
  poi_y <- poi$y

  if (fit_type == "linear") {
    m <- params[1]
    y_pred <- m * (x - poi_x) + poi_y
    constr_cost <- 0
  } else {  # Parabolic
    a <- params[1]
    b <- params[2]
    c <- params[3]
    constraints <- nl_constraints(params, poi_x, poi_y, x)
    constr_cost <- 100 * constraints$pass_poi + 100 * constraints$pos_grad
    y_pred <- a * x^2 + b * x + c
  }

  residuals <- y - y_pred
  l2_cost <- sum(residuals^2)
  return(l2_cost + constr_cost)  # Value to minimize
}

# Define a function to fit a profile to the left or right of a point of interest
fit_profile <- function(x, y, poi, slope_initial, fit_type = "linear") {
  initial_params <- if (fit_type == "linear") {
    c(slope_initial)
  } else {
    c(a = 0.5, b = 2, c = poi$y)
  }

  optim_result <- stats::optim(
    par = initial_params,
    fn = objective_function,
    x = x,
    y = y,
    poi = poi,
    fit_type = fit_type,
    method = "L-BFGS-B"
  )

  return(list(residual = optim_result$value, params = optim_result$par))
}

# Define a function to fit two splines around a point of interest
fit <- function(data, poi, fit_type = "linear") {
  x <- data$x
  y <- data$y

  poi_x <- poi$x
  poi_y <- poi$y

  # Left fit
  left_indices <- which(x <= poi_x)
  slope_initial_left <- (poi_y - y[left_indices][1]) / (poi_x - x[left_indices][1])
  result_left <- fit_profile(x = x[left_indices], y = y[left_indices], poi = poi, slope_initial = slope_initial_left, fit_type = fit_type)

  # Right fit
  right_indices <- which(x > poi_x)
  slope_initial_right <- (y[right_indices][length(right_indices)] - poi_y) / (x[right_indices][length(right_indices)] - poi_x)
  result_right <- fit_profile(x = x[right_indices], y = y[right_indices], poi = poi, slope_initial = slope_initial_right, fit_type = fit_type)

  total_residuals <- result_left$residual + result_right$residual

  return(list(residual = total_residuals, left_params = result_left$params, right_params = result_right$params))
}

# Define the function to seek the point of inflection
seek_inflection <- function(data, roi, step_x = 0.05, step_y = 0.1, fit_type = "linear") {
  grid_points <- make_grid(roi, step_x, step_y)
  best_residual <- Inf
  best_point <- NULL
  best_params_left <- NULL
  best_params_right <- NULL

  for (i in seq_len(nrow(grid_points))) {
    poi <- grid_points[i, ]
    result <- fit(data, poi, fit_type)
    if (result$residual < best_residual) {
      best_residual <- result$residual
      best_point <- poi
      best_params_left <- result$left_params
      best_params_right <- result$right_params
    }
  }

  return(list(inflection_point = best_point, left_params = best_params_left, right_params = best_params_right))
}
