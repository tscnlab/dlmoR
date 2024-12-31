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
  poi_diff <- a * poi_x**2 + b * poi_x + c - poi_y

  # Constraint 2: Ensure positive slope (dy/dx > 0) for all points in ascending segment
  slope_values <- 2 * a * x + b
  slope_min <- min(slope_values)

  if (slope_min >= 0) {
    slope_min <- 0
  }
  return(list(pass_poi = poi_diff**2, pos_grad = slope_min**2))
}

nl_constraints_new <- function(params, poi_x, poi_y, x) {
  a <- params[1]
  b <- params[2]
  c <- params[3]

  # Constraint 1: Parabola must pass through the POI
  poi_diff <- a * poi_x**2 + b * poi_x + c - poi_y

  # Constraint 2: Ensure positive slope (dy/dx > 0) for all points in ascending segment
  slope_values <- 2 * a * x + b
  slope_min <- tail(slope_values, n=1) # CHANGE (I dont think this had an impact). We only need to constraint the last point. It is gauranteed if this slope is >= and the one at the inflection point (other constraint), everything in between is.

  if (slope_min >= 0) {
    slope_min <- 0
  }
  return(list(pass_poi = poi_diff**2, pos_grad = slope_min**2))
}

poi_constraints<-function(params, slope_bounds, poi_x, fit_type = "linear"){
  if(fit_type == "linear"){
    slope_value<- params[1]
  }
  else{
    a <- params[1]
    b <- params[2]
    slope_value <- 2 * a * poi_x + b
  }

  #Constraint to ensure that slope at poi meets all conditions
    if (dplyr::between(slope_value, slope_bounds[1],slope_bounds[2])){
    slope_value <- 0
    } else{
      slope_value <- min(c(slope_bounds[1]-slope_value)**2,(slope_bounds[2]-slope_value)**2)
    }
    return(slope_value)
}

# Define objective function for constraining base fit edges along y axis
base_constraint <- function(y, min_y, poi_y, threshold = 2.3){ #TODO bring in threshold externally, not hardcoded
  if (dplyr::between(y[1],0,threshold)){
  left_constr<-0
}else{
  left_constr<-min((y[1]-threshold)**2,y[1]**2)
}
  if (dplyr::between(poi_y,min_y,threshold)){
    right_constr<-0
  }else{
    right_constr<-min((poi_y-threshold)**2,(poi_y-min_y)**2)
  }
  return(left_constr+right_constr)
}
# Define the objective function for line or parabola fitting
objective_function <- function(params, x, y, poi, fit_type, slope_bounds, base_id, region, threshold = threshold) {
  poi_x <- poi$x
  poi_y <- poi$y

  if (fit_type == "linear" & region =="ascending") {
    m <- params[1]
    y_pred <- m * (x - poi_x) + poi_y
    constr_cost <- 0
  } else if (fit_type == "linear"){ #base
    m <- params[1]
    y_pred <- m * (x - poi_x) + poi_y
    constr_cost <- base_constraint(y_pred, min(y), poi_y, threshold = threshold)*100
  }
  else {  # Parabolic
    a <- params[1]
    b <- params[2]
    c <- params[3]
    constraints <- nl_constraints(params, poi_x, poi_y, x)
    constr_cost <- 100 * constraints$pass_poi + 100 * constraints$pos_grad
    y_pred <- a * x**2 + b * x + c
  }
  poi_cost <- poi_constraints(params, slope_bounds, poi_x, fit_type)
  constr_cost <-constr_cost + poi_cost*100
  residuals <- y - y_pred
  residuals <- residuals*(1-base_id)+0.5*residuals*base_id #give base points 1/2 influence on fit
  residuals <- 0.5*residuals # give points on left 1/2 influence on fit as points on right
  l2_cost <- mean(residuals**2) # considering switching to mean so that fit is agnostic of number of data points
  return(l2_cost + constr_cost)  # Value to minimize
}

# Define the objective function for line or parabola fitting
# CHANGE New objective function for new fitting
objective_function_new <- function(params, x, y, poi, fit_type, slope_bounds, base_id, region, threshold = threshold) {
  poi_x <- poi$x
  poi_y <- poi$y

  if (fit_type == "linear" & region =="ascending") {
    m <- params[1]
    y_pred <- m * (x - poi_x) + poi_y
    constr_cost <- 0
  } else if (fit_type == "linear"){ #base
    m <- params[1]
    y_pred <- m * (x - poi_x) + poi_y
    constr_cost <- base_constraint(y_pred, min(y), poi_y, threshold = threshold)*100
  }
  else {  # Parabolic
    a <- params[1]
    b <- params[2]
    c <- poi$y  - a * poi$x**2 - b * poi$x # CHANGE we can actually compute C
    params <- c(a,b,c)
    constraints <- nl_constraints_new(params, poi_x, poi_y, x)
    constr_cost <- 100 * constraints$pos_grad # CHANGE now we dont need to constraint the optimization to go through the inflection point
    y_pred <- a * x**2 + b * x + c
  }
  poi_cost <- poi_constraints(params, slope_bounds, poi_x, fit_type)
  constr_cost <-constr_cost + poi_cost*100
  residuals <- y - y_pred
  residuals <- residuals*(1-base_id)+0.5*residuals*base_id #give base points 1/2 influence on fit # nolint: line_length_linter.
  #residuals <- 0.5*residuals # give points on left 1/2 influence on fit as points on right
  l2_cost <- mean(residuals**2) # considering switching to mean so that fit is agnostic of number of data points
  return(l2_cost + constr_cost)  # Value to minimize
}

# Define a function to fit a profile to the base or ascending of a point of interest
fit_profile <- function(x, y, poi, slope_initial, fit_type = "linear", region = "base", base_tangent = NULL, ascending_linear_tangent = NULL, base_id = NULL, threshold = threshold) {
  # initial_params <- if (fit_type == "linear") {
  #   c(slope_initial)
  #   #c(0)
  # } else {
  #   # c(a = 0.5, b = 2, c = poi$y)
  #   c(a = 0, b = ascending_linear_tangent, c = (ascending_linear_tangent*-poi$x)+poi$y)
  #   # c(a = -3, b = 160, c = -2000)
  # }

  if(region == "base"){
    slope_bounds <- c(-.2,.2)
    edge_bounds <- list(left = c(0,threshold), right = c(min(y), threshold))
    return(fit_linear(x,y,poi,base_id, slope_bounds, edge_bounds))
  } else if(region == "ascending" & fit_type == "linear") { #ascending
    slope_bounds <- c(-Inf, Inf)
    return(fit_linear(x,y,poi,base_id, slope_bounds))
  } else{ # parabolic fit on right side
    slope_lowerbound <- max(c(0,base_tangent,0.5*ascending_linear_tangent))
    slope_upperbound <- Inf
    slope_bounds <- c(slope_lowerbound, slope_upperbound)

    # CHANGE improved parabola fitting with only two paramters (requires also less constraints.)
    # I couldnt plot it, but the parameters and residuals looked more accurate when printed.
    initial_params<- c(a = 0, b = ascending_linear_tangent)
    optim_result <- stats::optim(
      par = initial_params,
      fn = objective_function_new,
      x = x,
      y = y,
      poi = poi,
      fit_type = fit_type,
      slope_bounds = slope_bounds,
      base_id = base_id,
      region = region,
      threshold = threshold,
      method = "L-BFGS-B"
    )
    a <- optim_result$par[1]
    b <- optim_result$par[2]
    c <- poi$y  - a * poi$x**2 - b * poi$x

    # CHANGE only care about the actual difference in prediction and truth. Dont care about any constraints (for search).

    y_pred <- a * x ** 2 + b * x + c
    delta <- y_pred - y
    residual <- delta ** 2
    residual <- residual*(1-base_id)+0.25*residual*base_id
    return(list(residual = sum(residual), params = list(a=a, b=b,c=c)))
  }

}

# perform linear fit
fit_linear <- function(x, y, poi, base_id = NULL, slope_bounds = NULL, edge_bounds = NULL){
  x_diff<- x - poi$x
  y_diff<- y - poi$y
  w<- (1-base_id)+0.5*base_id
  num<- sum(w*x_diff*y_diff)
  denom<- sum(w*x_diff**2)
  slope<- num/denom
  if(slope<min(slope_bounds)){
    slope<-min(slope_bounds)
  }else if(slope>max(slope_bounds)){
    slope<-max(slope_bounds)
  }
  if(!is.null(edge_bounds)){ #left side fit
    left_edge<-slope*(x[1]-poi$x)+poi$y
    right_edge<-poi$y

    # left edge bounds
    if(!dplyr::between(left_edge,edge_bounds$left[1],edge_bounds$left[2])){
      if(left_edge<edge_bounds$left[1]){
        left_edge<-edge_bounds$left[1]
      }else{
        left_edge<-edge_bounds$left[2]
      }
      slope<-(left_edge-right_edge)/(x[1]-poi$x)
    }

    # right edge bounds
    if(!dplyr::between(right_edge,edge_bounds$right[1],edge_bounds$right[2])){
      if(right_edge<edge_bounds$right[1]){
        right_edge<-edge_bounds$right[1]
      }else{
        right_edge<-edge_bounds$right[2]
      }
      slope<-(left_edge-right_edge)/(x[1]-poi$x)
    }
  }
  y_pred<-slope*(x-poi$x)+poi$y
  residuals <- y-y_pred
  residuals <- residuals*(1-base_id)+0.5*residuals*base_id #give base points 1/2 influence on fit
  # list(residual = optim_result$value, params = optim_result$par)
  return(list(residual = mean(residuals**2), params = c(slope)))
}

# Define a function to fit two splines around a point of interest
fit <- function(data, poi, fit_type = "linear", threshold = threshold) {
  x <- posixct_to_decimal(data$datetime, data$datetime)
  y <- data$melatonin
  base_id<-data$base
  #print("base_id")
  #print(base_id)

  poi_x <- poi$x
  poi_y <- poi$y

  # Base fit
  base_indices <- which(x <= poi_x)
  # slope_initial_base <- (poi_y - y[base_indices][1]) / (poi_x - x[base_indices][1])
  slope_initial_base <- 0

  result_base <- fit_profile(x = x[base_indices], y = y[base_indices], poi = poi, slope_initial = slope_initial_base, fit_type = "linear", region = "base", base_id = base_id[base_indices], threshold = threshold)

  # Ascending fit
  ascending_indices <- which(x > poi_x)
  # print("x")
  # print(x)
  # print("ascending")
  # print(ascending_indices)
  # print("poi_X")
  # print(poi_x)

  slope_initial_ascending <- (y[ascending_indices][length(ascending_indices)] - poi_y) / (x[ascending_indices][length(ascending_indices)] - poi_x)
  result_ascending <- fit_profile(x = x[ascending_indices], y = y[ascending_indices], poi = poi, slope_initial = slope_initial_ascending, fit_type = "linear", region = "ascending", base_id = base_id[ascending_indices], threshold = threshold)
  # print("initial ascending linear fit")
  # print(result_ascending)

  ascending_indices <- which(x > poi_x)
  if(length(ascending_indices) > 2){
  slope_initial_ascending<- result_ascending$params[1]
  # slope_initial_ascending<-(y[ascending_indices][length(ascending_indices)] - poi_y) / (x[ascending_indices][length(ascending_indices)] - poi_x)
  result_ascending <- fit_profile(x = x[ascending_indices], y = y[ascending_indices], poi = poi, slope_initial = slope_initial_ascending, fit_type = "parabolic", base_tangent = -result_base$params[1], ascending_linear_tangent = result_ascending$params[1], region = "ascending", base_id = base_id[ascending_indices], threshold = threshold)
  }
  total_residuals <- result_base$residual + result_ascending$residual
  # total_residuals <- result_ascending$residual
  # print("result_base")
  # print(result_base)
  # print("result_ascending")
  # print(result_ascending)
  return(list(residual = total_residuals, base_params = result_base$params, ascending_params = result_ascending$params))
}

# Define the function to seek the point of inflection
#seek_inflection <- function(data, roi, step_x = 0.05, step_y = 0.1, fit_type = "linear") {
seek_inflection <- function(data, threshold = threshold, roi, step_x = 0.1, step_y = 0.2, step_size_small = 0.01, fit_type = "linear") {

  grid_points <- make_grid(roi, step_x, step_y)
  best_residual <- Inf
  best_point <- NULL
  best_params_base <- NULL
  best_params_ascending <- NULL
  res_big<-NULL
  grid_big<-grid_points

  for (i in seq_len(nrow(grid_points))) {
    poi <- grid_points[i, ]
    result <- fit(data, poi, fit_type, threshold = threshold)
    res_big[i]<-result$residual
    # print("res")
    # print(res)
    if (result$residual < best_residual) {
      best_residual <- result$residual
      best_point <- poi
      best_params_base <- result$base_params
      best_params_ascending <- result$ascending_params
    }
  }


  # CHANGE: Search of an area that includes the 10% smallest points (probably what the paper wants)
  #best_10_per <- order(res_big)[1:as.integer(0.1*nrow((grid_points)))]
  #best_points <- grid_points[best_10_per,]
  #print(res_big[best_10_per][1:10])
  #print(best_points[1:10,])
  #min_y <- max(min(min(best_points$y), 2.3), 0)
  #max_y <- max(min(max(best_points$y) + step_y, 2.3), 0)
  #min_x <- min(best_points$x) - step_x
  #max_x <- max(best_points$x) + step_x

  # CHANGE search area around best one. Gives you the result you want though
  # unsure about the min melatonin. That was a constraint right.
  # best_10_per <- order(res_big)[1:10]
  # best_points <- grid_points[best_10_per, ][1:1,]
  # print(res_big[best_10_per][1:6])
  # print(best_points[1:6,])
  # print(best_point[1,])
  # min_y <- max(min(min(best_points$y) - step_y, 2.3), min(data$melatonin))
  # max_y <- max(min(max(best_points$y) + step_y, 2.3), min(data$melatonin))
  # min_x <- min(best_points$x) - step_x
  # max_x <- max(best_points$x) + step_x

  # print(min_x)
  # print(min_y)
  # print(max_x)
  # print(max_y)
  # roi <- list(x = c(min_x, max_x), y = c(min_y, max_y))
  # grid_points <- make_grid(roi, 0.01, 0.01)
  # grid_small <- grid_points
  grid_points <- reduce_grid2(res_big = res_big, grid_big = grid_big, step_x = step_x, step_y = step_y, step_size_small = step_size_small, threshold = threshold)
  grid_small <- grid_points
  res_small<-NULL

  # for (i in seq_len(nrow(grid_points))) {
  #   poi <- grid_points[i, ]
  #   result <- fit(data, poi, fit_type, threshold = threshold)
  #   res_small[i]<-result$residual
  #   if (result$residual < best_residual) {
  #     best_residual <- result$residual
  #     best_point <- poi
  #     best_params_base <- result$base_params
  #     best_params_ascending <- result$ascending_params
  #   }
  # }

  # print(sort(res_small)[1:10])
  # best_10_per <- order(res_small)[1:10]
  # best_points <- grid_points[best_10_per,]
  # print(best_points)
  return(list(inflection_point = best_point, base_params = best_params_base, ascending_params = best_params_ascending, grid_big = grid_big, res_big = res_big, grid_small = grid_small, res_small = res_small))
}

reduce_grid1 <- function(res_big,grid_big, step_size_small, threshold = threshold){
  best_10_per <- order(res_big)[1:10]
  # best_points <- grid_points[best_10_per, ][1:1,]
  best_points <- grid_big[best_10_per, ]
  print(best_points)
  # min_y <- max(min(min(best_points$y) - 0.2, 2.3), min(data$melatonin)) #TODO threshold not 2.3
  # max_y <- max(min(max(best_points$y) + 0.2, 2.3), min(data$melatonin))
  min_y <- min(min(best_points$y) - 0.2, 0) #TODO threshold not 2.3
  max_y <- min(max(best_points$y) + 0.2, threshold)
  min_x <- min(best_points$x) - 0.1
  max_x <- max(best_points$x) + 0.1
  print(min_y)
  print(max_y)
  print(min_x)
  print(max_x)
  roi <- list(x = c(min_x, max_x), y = c(min_y, max_y))
  grid_points <- make_grid(roi, step_size_small, step_size_small)
}

reduce_grid2 <- function(res_big, grid_big, step_x, step_y, step_size_small, threshold = threshold){
  best_10_per <- order(res_big)[1:as.integer(0.1*nrow((grid_big)))]
  best_points <- grid_big[best_10_per,]
  print(best_points)
  #print(res_big[best_10_per][1:10])
  #print(best_points[1:10,])
  # min_y <- max(min(min(best_points$y) - 0.2, 2.3), 0)
  # max_y <- max(min(max(best_points$y) + 0.2, 2.3), 0)
  # min_x <- min(best_points$x) - 0.1
  # # max_x <- max(best_points$x) + 0.1
  # min_y <- min(best_points$y) - step_y
  # max_y <- max(best_points$y) + step_y #TODO threshold not 2.3
  # min_x <- min(best_points$x) - step_x
  # max_x <- max(best_points$x) + step_x

  min_y <- max(min(best_points$y) - step_y, min(grid_big$y))
  max_y <- min(max(best_points$y) + step_y, max(grid_big$y)) #TODO threshold not 2.3
  min_x <- min(best_points$x) - step_x
  max_x <- max(best_points$x) + step_x
  print(min_y)
  print(max_y)
  print(min_x)
  print(max_x)

  roi <- list(x = c(min_x, max_x), y = c(min_y, max_y))
  print(roi)
  grid_points <- make_grid(roi, step_size_small, step_size_small)
}

# run this script to get inflection
get_inflection <- function(profile_data, threshold = 2.3, posix_roi, fit_type = "linear"){
  roi<-list(x = posixct_to_decimal(c(posix_roi$x_start, posix_roi$x_end), profile_data$datetime), y = c(posix_roi$y_min, posix_roi$y_max))
  if ("intermediate"%in%colnames(profile_data)){
  poi<-seek_inflection(dplyr::filter(profile_data,base == 1 | ascending == 1 | intermediate == 1), threshold = threshold, roi, fit_type = fit_type)
  }else{
    poi<-seek_inflection(dplyr::filter(profile_data,base == 1 | ascending == 1), threshold = threshold, roi, fit_type = fit_type)
  }
  return(poi)
}


# load('/Users/langert1/Library/CloudStorage/OneDrive-AaltoUniversity/Documents/DLMO/20decenviro.RData')
# source('time_to_decimal.R')
#ip<-get_inflection(dlmo204FDd2v4$prof, dlmo204FDd2v4$roi)
#print(ip$inflection_point)
