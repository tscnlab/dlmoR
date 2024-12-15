# Load required libraries
#if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
#if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

# Function to convert POSIXct times to decimal hours
# posixct_to_decimal <- function(posix_times, profile_datetime) {
#   posix_times <- as.POSIXct(posix_times, tz = "UTC")
#   posix_origin <- profile_datetime[1]
#   origin_date <- as.Date(posix_origin)
#
#   days_elapsed <- as.numeric(as.Date(posix_times) - origin_date)
#   hours <- as.numeric(format(posix_times, "%H"))
#   minutes <- as.numeric(format(posix_times, "%M"))
#
#   decimal_time_today <- hours + (minutes / 60)
#   decimal_hours <- (days_elapsed * 24) + decimal_time_today
#   return(decimal_hours)
# }

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

  # Constraint 2: Ensure positive slope (dy/dx > 0) for all points in ascending segment
  slope_values <- 2 * a * x + b
  slope_min <- min(slope_values)

  if (slope_min >= 0) {
    slope_min <- 0
  }

  return(list(pass_poi = poi_diff^2, pos_grad = slope_min^2))
}

poi_constraints<-function(params, slope_bounds, poi_x, fit_type = "linear"){
  if(fit_type == "linear"){
    slope_value<- params[1]
  }
  else{
    a <- params[1]
    b <- params[2]
    slope_value <- 2 * a * poi_x + b
    #slope_value <- 2 * a * 0 + b
  }


  #Constraint to ensure that slope at poi meets all conditions
    if (dplyr::between(slope_value, slope_bounds[1],slope_bounds[2])){
    slope_value <- 0
    } else{
      slope_value <- min(c(slope_bounds[1]-slope_value)^2,(slope_bounds[2]-slope_value)^2)
    }
    return(slope_value)
}

# Define objective function for constraining base fit edges along y axis
base_constraint <- function(y, min_y, poi_y, threshold = 2.3){ #TODO bring in threshold externally, not hardcoded
  if (dplyr::between(y[1],0,threshold)){
  left_constr<-0
}else{
  left_constr<-min((y[1]-threshold)^2,y[1]^2)
}
  if (dplyr::between(poi_y,min_y,threshold)){
    right_constr<-0
  }else{
    right_constr<-min((poi_y-threshold)^2,(poi_y-min_y)^2)
  }
  return(left_constr+right_constr)
}
# Define the objective function for line or parabola fitting
objective_function <- function(params, x, y, poi, fit_type, slope_bounds, base_id, region) {
  poi_x <- poi$x
  poi_y <- poi$y

  if (fit_type == "linear" & region =="ascending") {
    m <- params[1]
    y_pred <- m * (x - poi_x) + poi_y
    constr_cost <- 0
  } else if (fit_type == "linear"){ #base
    m <- params[1]
    y_pred <- m * (x - poi_x) + poi_y
    constr_cost <- base_constraint(y_pred, min(y), poi_y)*100
  }
  else {  # Parabolic
    a <- params[1]
    b <- params[2]
    c <- params[3]
    constraints <- nl_constraints(params, poi_x, poi_y, x)
    constr_cost <- 100 * constraints$pass_poi + 100 * constraints$pos_grad
    y_pred <- a * x^2 + b * x + c
  }
  poi_cost <- poi_constraints(params, slope_bounds, poi_x, fit_type)
  constr_cost <-constr_cost + poi_cost*100
  residuals <- y - y_pred
  # print(length(residuals))
  # print(length(base_id))
  residuals <- residuals*(1-base_id)+0.5*residuals*base_id #give base points 1/2 influence on fit
  l2_cost <- mean(residuals^2) # considering switching to mean so that fit is agnostic of number of data points
  return(l2_cost + constr_cost)  # Value to minimize
}

# Define a function to fit a profile to the base or ascending of a point of interest
fit_profile <- function(x, y, poi, slope_initial, fit_type = "linear", region = "base", base_tangent = NULL, ascending_linear_tangent = NULL, base_id = NULL) {
  initial_params <- if (fit_type == "linear") {
    c(slope_initial)
    #c(0)
  } else {
    # c(a = 0.5, b = 2, c = poi$y)
    c(a = 0, b = ascending_linear_tangent, c = (ascending_linear_tangent*-poi$x)+poi$y)
    # c(a = -3, b = 160, c = -2000)
  }

  if(region == "base"){
    slope_bounds = c(-.2,.2)
  } else if(region == "ascending" & fit_type == "linear") { #ascending
    slope_bounds = c(-Inf, Inf)
  } else{
    slope_lowerbound = max(c(0,base_tangent,0.5*ascending_linear_tangent))
    # print("base tangent")
    # print(base_tangent)
    # print("1/2 ascending")
    # print(0.5*ascending_linear_tangent)
    # print("slope lower bound")
    # print(slope_lowerbound)
    slope_upperbound = Inf
    slope_bounds = c(slope_lowerbound, slope_upperbound)
  }
  optim_result <- stats::optim(
    par = initial_params,
    fn = objective_function,
    x = x,
    y = y,
    poi = poi,
    fit_type = fit_type,
    slope_bounds = slope_bounds,
    base_id = base_id,
    region = region,
    method = "L-BFGS-B"
  )
  return(list(residual = optim_result$value, params = optim_result$par))
}


# Define a function to fit two splines around a point of interest
fit <- function(data, poi, fit_type = "linear") {
  x <- posixct_to_decimal(data$datetime, data$datetime)
  y <- data$melatonin
  base_id<-data$base
  poi_x <- poi$x
  poi_y <- poi$y

  # Base fit
  base_indices <- which(x <= poi_x)
  # slope_initial_base <- (poi_y - y[base_indices][1]) / (poi_x - x[base_indices][1])
  slope_initial_base <- 0
  result_base <- fit_profile(x = x[base_indices], y = y[base_indices], poi = poi, slope_initial = slope_initial_base, fit_type = "linear", region = "base", base_id = base_id[base_indices])

  # Ascending fit
  ascending_indices <- which(x > poi_x)
  slope_initial_ascending <- (y[ascending_indices][length(ascending_indices)] - poi_y) / (x[ascending_indices][length(ascending_indices)] - poi_x)
  result_ascending <- fit_profile(x = x[ascending_indices], y = y[ascending_indices], poi = poi, slope_initial = slope_initial_ascending, fit_type = "linear", region = "ascending", base_id = base_id[ascending_indices])
  # print("initial ascending linear fit")
  # print(result_ascending)

  ascending_indices <- which(x > poi_x)
  if(length(ascending_indices) > 2){
  slope_initial_ascending<- result_ascending$params[1]
  # slope_initial_ascending<-(y[ascending_indices][length(ascending_indices)] - poi_y) / (x[ascending_indices][length(ascending_indices)] - poi_x)
  result_ascending <- fit_profile(x = x[ascending_indices], y = y[ascending_indices], poi = poi, slope_initial = slope_initial_ascending, fit_type = "parabolic", base_tangent = -result_base$params[1], ascending_linear_tangent = result_ascending$params[1], region = "ascending", base_id = base_id[ascending_indices])
  }
  total_residuals <- result_base$residual + result_ascending$residual
  # total_residuals <- result_ascending$residual
  return(list(residual = total_residuals, base_params = result_base$params, ascending_params = result_ascending$params))
}

# Define the function to seek the point of inflection
seek_inflection <- function(data, roi, step_x = 0.05, step_y = 0.1, fit_type = "linear") {
  grid_points <- make_grid(roi, step_x, step_y)
  best_residual <- Inf
  best_point <- NULL
  best_params_base <- NULL
  best_params_ascending <- NULL
  res<-NULL

  for (i in seq_len(nrow(grid_points))) {
    poi <- grid_points[i, ]
    result <- fit(data, poi, fit_type)
    # res[i]<-result$residual
    # print("res")
    # print(res)
    if (result$residual < best_residual) {
      best_residual <- result$residual
      best_point <- poi
      best_params_base <- result$base_params
      best_params_ascending <- result$ascending_params
    }
  }
  print("best_residual")
  print(best_residual)
  return(list(inflection_point = best_point, base_params = best_params_base, ascending_params = best_params_ascending))
}

# run this script to get inflection
get_inflection <- function(profile_data, posix_roi, fit_type = "linear"){
  roi<-list(x = posixct_to_decimal(c(posix_roi$x_start, posix_roi$x_end), profile_data$datetime), y = c(posix_roi$y_min, posix_roi$y_max))
  poi<-seek_inflection(dplyr::filter(profile_data,base ==1 | ascending == 1), roi, fit_type = fit_type)
  return(poi)
}
