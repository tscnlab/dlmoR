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



#TL: REMOVED OBJECTIVE_FUNCTION WE ONLY USE _NEW, SO PREVENTING CONSFUSING.

# Define the objective function for line or parabola fitting
# CHANGE New objective function for new fitting
objective_function_new <- function(params, x, y, poi, fit_type, slope_bounds, base_id, region, weight_base) {
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
    c <- poi$y  - a * poi$x**2 - b * poi$x # CHANGE we can actually compute C
    params <- c(a,b,c)
    constraints <- nl_constraints_new(params, poi_x, poi_y, x)
    constr_cost <- 100 * constraints$pos_grad # CHANGE now we dont need to constraint the optimization to go through the inflection point
    y_pred <- a * x**2 + b * x + c
  }
  poi_cost <- poi_constraints(params, slope_bounds, poi_x, fit_type)
  constr_cost <-constr_cost + poi_cost*100
  residuals <- y - y_pred
  # print(length(residuals))
  # print(length(base_id))

  if (weight_base){
  residuals <- residuals*(1-base_id)+0.5*residuals*base_id #TL: NO CHANGE BUT HIGHLIGHT. HERE WE WEIGHT TO MAKE SURE ASCENDING COUNTS MORE DURING FITTING OF PARABOLA. 
  }
  l2_cost <- mean(residuals**2) # considering switching to mean so that fit is agnostic of number of data points
  return(l2_cost + constr_cost)  # Value to minimize
}

# Define a function to fit a profile to the base or ascending of a point of interest
fit_profile <- function(x, y, poi, fit_type = "linear", region = "base", base_tangent = NULL, ascending_linear_tangent = NULL, base_id = NULL, weight_base=TRUE) {
  #TL: REMOVE IF7ELSE THAT SET INITIALIZATIONS. WE ARE LITERALLY NOT USING IT. 
  if(region == "base"){
    slope_bounds <- c(-.2,.2)
    edge_bounds <- list(left = c(0,2.3), right = c(min(y), 2.3))
    result <- fit_linear(x,y,poi,base_id, slope_bounds, edge_bounds, weight_base = weight_base)
    params <- result$params
    y_pred<-params[1]*(x-poi$x)+poi$y
    #TL: REMOVED RESIDUAL, RECOMPUATING BELOW. 

  } else if(region == "ascending" & fit_type == "linear") { #ascending
    slope_bounds <- c(-Inf, Inf)
    result <- fit_linear(x,y,poi,base_id, slope_bounds, weight_base = weight_base)
    params <- result$params
    y_pred <- params[1]*(x-poi$x)+poi$y 
    # TL: REMOVED RETURN, RECOMPUTING RESIDUAL BELOW. 

  } else{ # parabolic fit on right side
    slope_lowerbound <- max(c(0,base_tangent,0.5*ascending_linear_tangent))
    slope_upperbound <- Inf
    slope_bounds <- c(slope_lowerbound, slope_upperbound)

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
      weight_base = weight_base,
      method = "L-BFGS-B"
    )
    a <- optim_result$par[1]
    b <- optim_result$par[2]
    c <- poi$y  - a * poi$x**2 - b * poi$x
    y_pred <- a * x ** 2 + b * x + c
    params <- list(a=a, b=b,c=c)
    #TL: TOOK OUT THE RESIDUAL COMPUTATION OF THE IF-ELSE STATEMENT SO WE CAN USE IT FOR ALL 3 SCENARIOS IN THIS IF, ELSE IF, ELSE STATEMENT
  }

  delta <- y_pred - y
  print('---')
  print(region)
  print(fit_type)
  print(base_id)
  print(delta)
  if(weight_base){
  residual <- delta*(1-base_id)+0.5*delta*base_id #TL: HERE WE WEIGHT. THIS IS THE RESIDUAL THAT GETS RETURN TO FIT AND VIA THAT TO INFLECTION POINT SEARCH
  }
  else{
  residual <- delta
  }
  print(residual)
  return(list(residual = mean(residual ** 2), params=params))
}
# perform linear fit
fit_linear <- function(x, y, poi, base_id = NULL, slope_bounds = NULL, edge_bounds = NULL, weight_base=TRUE){
  x_diff<- x - poi$x
  y_diff<- y - poi$y
  if(weight_base) {
  w <- (1-base_id)+0.5*base_id
  }
  else {
  w <- 1
  }
  num<- sum(w*x_diff*y_diff) #TL: NO CHANGE BUT HIGHLIGHT HERE WE WEIGHT to make sure that our linear fit takes ascending points more into account.
  denom<- sum(w*x_diff**2) #TL: NO CHANGE BUT HIGHLIGHT. HERE WE WEIGHT to make sure that our linear fit takes ascending points more into account.
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

  ## THERE WAS SOME RESIDUAL COMPUTATOIN HERE. I REMOVED IT BECAUSE IN THE NEW VERSION IT IS NOT NEEDED, WANTED TO PREVENT CONFUSION
  return(list(params = c(slope)))
}

# Define a function to fit two splines around a poinÍt of interest
fit <- function(data, poi, fit_type = "linear") {
  x <- posixct_to_decimal(data$datetime, data$datetime)
  y <- data$melatonin
  base_id<-data$base

  poi_x <- poi$x
  poi_y <- poi$y
  weight_base <- TRUE #TL TURN OFF WEIGHTING HERE
  # Base fit
  left_indcs <- which(x <= poi_x)
  result_base <- fit_profile(x = x[left_indcs], y = y[left_indcs], poi = poi, fit_type = "linear", region = "base", base_id = base_id[left_indcs], weight_base=weight_base)

  # Ascending fit
  right_indcs <- which(x > poi_x)
  result_ascending <- fit_profile(x = x[right_indcs], y = y[right_indcs], poi = poi, fit_type = "linear", region = "ascending", base_id = base_id[right_indcs], weight_base=weight_base)

  if(length(right_indcs) > 2){
    slope_initial_ascending <- result_ascending$params[1]
    result_ascending <- fit_profile(x = x[right_indcs], y = y[right_indcs], poi = poi, fit_type = "parabolic", base_tangent = -result_base$params[1], ascending_linear_tangent = result_ascending$params[1], region = "ascending", base_id = base_id[right_indcs], weight_base=weight_base)
  }
  total_residuals <- result_base$residual + result_ascending$residual
  return(list(residual = total_residuals, base_params = result_base$params, ascending_params = result_ascending$params))
}

# Define the function to seek the point of inflection
seek_inflection <- function(data, roi, step_x = 0.05, step_y = 0.1, fit_type = "linear") {
  grid_points <- make_grid(roi, step_x, step_y)
  best_residual <- Inf
  best_point <- NULL
  best_params_base <- NULL
  best_params_ascending <- NULL
  res_big<-NULL
  grid_big<-grid_points

  for (i in seq_len(nrow(grid_points))) {
    poi <- grid_points[i, ]
    result <- fit(data, poi, fit_type)
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


  best_10_per <- order(res_big)[1:10]
  best_points <- grid_points[best_10_per, ][1:1,]
  min_y <- max(min(min(best_points$y) - step_y, 2.3), min(data$melatonin))
  max_y <- max(min(max(best_points$y) + step_y, 2.3), min(data$melatonin))
  min_x <- min(best_points$x) - step_x
  max_x <- max(best_points$x) + step_x

  roi <- list(x = c(min_x, max_x), y = c(min_y, max_y))
  grid_points <- make_grid(roi, 0.01, 0.01)
  grid_small <- grid_points
  res_small<-NULL

  for (i in seq_len(nrow(grid_points))) {
    poi <- grid_points[i, ]
    result <- fit(data, poi, fit_type)
    res_small[i]<-result$residual
    if (result$residual < best_residual) {
      best_residual <- result$residual
      best_point <- poi
      best_params_base <- result$base_params
      best_params_ascending <- result$ascending_params
    }
  }

  print(sort(res_small)[1:10])
  print(best_residual)
  best_10_per <- order(res_small)[1:10]
  best_points <- grid_points[best_10_per,]
  return(list(inflection_point = best_point, base_params = best_params_base, ascending_params = best_params_ascending, grid_big = grid_big, res_big = res_big, grid_small = grid_small, res_small = res_small))
}

# run this script to get inflection
get_inflection <- function(profile_data, posix_roi, fit_type = "linear"){
  roi<-list(x = posixct_to_decimal(c(posix_roi$x_start, posix_roi$x_end), profile_data$datetime), y = c(posix_roi$y_min, posix_roi$y_max))
  if ("intermediate"%in%colnames(profile_data)){
    poi<-seek_inflection(dplyr::filter(profile_data,base == 1 | ascending == 1 | intermediate == 1), roi, fit_type = fit_type)
  }else{
    poi<-seek_inflection(dplyr::filter(profile_data,base == 1 | ascending == 1), roi, fit_type = fit_type)
  }
  return(poi)
}


load('/Users/langert1/Library/CloudStorage/OneDrive-AaltoUniversity/Documents/DLMO/20decenviro.RData')
source('time_to_decimal.R')
ip<-get_inflection(dlmo204FDd2v4$prof, dlmo204FDd2v4$roi)
print(ip$inflection_point)