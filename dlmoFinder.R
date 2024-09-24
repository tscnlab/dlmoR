# *** Parallelogram Fitter ***

library(R6)
library(ggplot2)


# Define the ParallelogramOptimizer class
ParallelogramOptimizer <- R6Class(
  "ParallelogramOptimizer",
  public = list(
    x = NULL,
    y = NULL,
    y0 = NULL,
    y1 = NULL,
    params = NULL,
    
    initialize = function() {
      self$params <- NULL
    },
    
    get_corners = function(x0, y0, x1, y1, slope) {
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
    },
    
    constraints = function(params) {
      x0 <- params[1]
      x1 <- params[2]
      slope <- params[3]
      constraint_vals <- numeric()
      
      corners <- self$get_corners(x0, self$y0, x1, self$y1, slope)
      lower_left <- corners[[1]]
      lower_right <- corners[[2]]
      upper_right <- corners[[3]]
      upper_left <- corners[[4]]
      
      for (i in seq_along(self$x)) {
        xi <- self$x[i]
        yi <- self$y[i]
        
        y_lower <- if (xi <= lower_right[1])
          self$y0
        else
          self$y0 + slope * (xi - lower_right[1])
        y_upper <- if (xi >= upper_left[1])
          self$y1
        else
          self$y1 + slope * (xi - upper_left[1])
        
        x_lower <- min(lower_left[1], upper_left[1])
        x_upper <- max(upper_right[1], lower_right[1])
        
        constraint_vals <- c(constraint_vals, y_upper - yi)
        constraint_vals <- c(constraint_vals, yi - y_lower)
        constraint_vals <- c(constraint_vals, xi - x_lower)
        constraint_vals <- c(constraint_vals, x_upper - xi)
      }
      return(constraint_vals)
    },
    
    objective = function(params) {
      x0 <- params[1]
      x1 <- params[2]
      slope <- params[3]
      corners <- self$get_corners(x0, self$y0, x1, self$y1, slope)
      v1 <- corners[[2]] - corners[[1]]
      v2 <- corners[[4]] - corners[[1]]
      area <- abs(v1[1] * v2[2] - v1[2] * v2[1])
      
      c_penalty <- min(self$constraints(params))
      if (c_penalty > 0) {
        c_penalty <- 0
      }
      return(area + 1e3 * c_penalty ^ 2)
    },
    
    optimize_parallelogram = function(x, y) {
      self$x <- x
      self$y <- y
      self$y0 <- min(y)
      self$y1 <- max(y)
      x0_initial <- min(self$x)
      x1_initial <- max(self$x)
      slope_initial <- (self$y1 - self$y0) / (x1_initial - x0_initial)
      
      initial_guess <- c(x0_initial, x1_initial, slope_initial)
      
      result <- optim(
        par = initial_guess,
        fn = self$objective,
        method = "L-BFGS-B",
        lower = c(-Inf, -Inf, 0),
        upper = c(Inf, Inf, Inf)
      )
      
      self$params <- result$par
      return(self$params)
    },
    
    plot_parallelogram = function(existing_plot = NULL,
                                  params = NULL) {
      # Use the provided params or fall back to self$params
      if (is.null(params)) {
        if (is.null(self$params)) {
          stop("Parameters are not optimized yet. Run parallelogram_truncation() first.")
        }
        params <- self$params
      }
      
      x0 <- params[1]
      x1 <- params[2]
      slope <- params[3]
      
      corners <- self$get_corners(x0, self$y0, x1, self$y1, slope)
      
      parallelogram_x <- c(corners[[1]][1], corners[[2]][1], corners[[3]][1], corners[[4]][1], corners[[1]][1])
      parallelogram_y <- c(corners[[1]][2], corners[[2]][2], corners[[3]][2], corners[[4]][2], corners[[1]][2])
      
      df <- data.frame(x = parallelogram_x, y = parallelogram_y)
      
      if (is.null(existing_plot)) {
        stop("No existing plot provided for overlay.")
      }
      
      # Overlay parallelogram on existing plot
      p <- existing_plot +
        geom_path(
          data = df,
          aes(x = x, y = y),
          color = "purple",
          size = 1
        ) +
        geom_segment(
          aes(
            x = corners[[1]][1],
            y = corners[[1]][2],
            xend = corners[[3]][1],
            yend = corners[[3]][2]
          ),
          color = "purple",
          linetype = "dashed"
        ) +
        geom_segment(
          aes(
            x = corners[[2]][1],
            y = corners[[2]][2],
            xend = corners[[4]][1],
            yend = corners[[4]][2]
          ),
          color = "purple",
          linetype = "dashed"
        )
      
      # p <- existing_plot
      
      print(p)
    }
  )
)

# **Region Fitter**:

# Load R6 package
library(R6)

# Define the InflectionPointSeeker class
RegionFitter <- R6Class(
  "RegionFitter",
  public = list(
    fit_type = NULL,
    
    # Constructor to initialize the object with parameters
    initialize = function() {
      
    },
    
    # Define the nonlinear constraints
    nl_constraints = function(params, poi_x, poi_y, x) {
      a <- params[1]
      b <- params[2]
      c <- params[3]
      
      # Constraint 1: Parabola must pass through the POI
      poi_diff <- a * poi_x ^ 2 + b * poi_x + c - poi_y
      
      # Constraint 2: Ensure positive slope (dy/dx > 0) for all points to the right of the POI
      slope_values <- 2 * a * x + b
      slope_min <- min(slope_values)
      
      if (slope_min >= 0) {
        slope_min <- 0
      }
      
      # Return the constraints: negative values indicate constraint violations
      return(list(pass_poi = poi_diff ^ 2, pos_grad = slope_min ^
                    2))
    },
    
    # Objective function to fit line with the constraint it pass through POI
    objective_function = function(params, x, y, poi, fit_type) {
      # Extract the slopes for the left and right lines
      
      # Extract poi coords
      poi_x <- poi[1, 1]
      poi_y <- poi[1, 2]
      
      if (fit_type == "linear") {
        # Calculate predicted values for fit lines
        m <- params[1]
        y_pred <- m * (x - poi_x) + poi_y
        constr_cost <- 0
      }
      else{
        #parabolic
        a <- params[1]
        b <- params[2]
        c <- params[3]
        constr <- self$nl_constraints(params, poi_x, poi_y, x)
        constr_cost <- 100 * constr$pass_poi + 100 *
          constr$pos_grad
        y_pred <- a * x ^ 2 + b * x + c
      }
      # Calculate residuals (difference between actual and predicted values)
      residuals <- y - y_pred
      
      # Sum of squared residuals
      l2_cost <- sum(residuals ^ 2)
      return(l2_cost + constr_cost)  # This is the value to minimize
    },
    
    # Method for a best fit line
    fit_profile = function(x, y, poi, slope_initial, fit_type =
                             "linear") {
      if (fit_type == "linear") {
        initial_params <- c(slope_initial)
      }
      else{
        initial_params <- c(a = 0.5, b = 2, c = poi[1, 2])
      }
      result <- optim(
        par = initial_params,
        fn = self$objective_function,
        x = x,
        y = y,
        poi = poi,
        fit_type = fit_type,
        method = "L-BFGS-B",
        # lower = c(-Inf),
        # upper = c(Inf)
      )
      
      # params <- result$par
      
      # residual <- self$objective_function(params, data_points = data , poi = poi)
      
      return (list(residual = result$value, para =
                     result$par))
    },
    
    # Method for fitting two splines
    fit = function(data, poi, fit_type = "linear") {
      x <- data$x
      y <- data$y
      
      poi_x <- poi[1, 1]
      poi_y <- poi[1, 2]
      
      # do left fit
      left_indices <- which(x <= poi_x)
      slope_initial <- (poi_y - y[left_indices][1]) / (poi_x - x[left_indices][1])
      
      result_left <- self$fit_profile(
        x = data$x[left_indices],
        y = data$y[left_indices],
        poi = poi,
        slope_initial = slope_initial
      )
      
      # do right fit
      right_indices <- which(x > poi_x)
      slope_initial <- (y[right_indices][length(right_indices)] - poi_y) /
        (x[right_indices][length(right_indices)] - poi_x)
      result_right <- self$fit_profile(
        x = data$x[right_indices],
        y = data$y[right_indices],
        poi = poi,
        slope_initial = slope_initial,
        fit_type = fit_type
      )
      
      # return summed cost
      total_residuals <- result_left$residual + result_right$residual
      return (
        list(
          residual = total_residuals,
          left_para = result_left$para,
          right_para = result_right$para
        )
      )
    }
  )
)

#TODO:
# - Add parabola
# - make base segments count half

# **Inflection Point Seeker**

# Load R6 package
library(R6)
library(ggplot2)

# Define the InflectionPointSeeker class
InflectionPointSeeker <- R6Class(
  "InflectionPointSeeker",
  public = list(
    roi = NULL,
    # Store the ROI object
    region_fitter = NULL,
    residuals = NULL,
    grid_points = NULL,
    
    # Constructor to initialize the object with parameters
    initialize = function() {
      self$region_fitter <- RegionFitter$new()
    },
    
    # Method to create the grid of points
    make_grid = function(roi, step_x, step_y) {
      # Extract xmin, xmax, ymin, ymax from ROI
      xmin <- roi$x[1]
      xmax <- roi$x[2]
      ymin <- roi$y[1]
      ymax <- roi$y[2]
      
      # Generate sequences for x and y coordinates
      x_seq <- seq(from = xmin, to = xmax, by = step_x)
      y_seq <- seq(from = ymin, to = ymax, by = step_y)
      
      # Create a grid of points
      grid_points <- expand.grid(x = x_seq, y = y_seq)
      
      # Return the grid points
      print('inflection point seeker: grid points')
      print(grid_points)
      return(grid_points)
    },
    
    seek = function(data,
                    roi,
                    step_x = 0.05,
                    step_y = 0.1,
                    fit_type = "linear") {
      self$grid_points <- self$make_grid(roi, step_x, step_y)
      best_residual <- Inf
      best_point <- NULL
      best_params_left <- NULL
      best_params_right <- NULL
      self$residuals <- c()
      for (i in 1:length(self$grid_points$x)) {
        poi <- self$grid_points[i, ]
        result <- self$region_fitter$fit(data =
                                           data,
                                         poi = poi,
                                         fit_type = fit_type)
        residual <- result$residual
        params_left <- result$left_para
        params_right <- result$right_para
        self$residuals[length(self$residuals) +
                         1] <- residual
        if (residual < best_residual) {
          best_residual <- residual
          best_point <- poi
          best_params_left <- params_left
          best_params_right <- params_right
          # print(best_params_left)
        }
      }
      return (
        list(
          inflection_point = best_point,
          left_para = best_params_left,
          right_para = best_params_right
        )
      )
    }
  )
)

# **DLMO Class**

library(R6)
library(ggplot2)

# Define the DLMO class
DLMO <- R6Class(
  "DLMO",
  public = list(
    input_data = NULL,
    threshold = NULL,
    trimmed_data = NULL,
    parallelogram_optimizer = ParallelogramOptimizer$new(),
    parallelogram_params = NULL,
    n_data = NULL,
    roi = NULL,
    coarse_roi = NULL,
    inflectionPointSeeker = InflectionPointSeeker$new(),
    result = NULL,
    
    # Constructor to initialize data and threshold
    initialize = function(input_data, threshold) {
      self$input_data <- input_data
      print('input_data')
      print(self$input_data)
      self$threshold <- threshold
      self$trimmed_data <- self$clean_data()
      print('trimmed_data')
      print(self$trimmed_data)
      self$n_data <- length(self$trimmed_data$x)  # Assuming 'x' represents the length of the data points
      
      # Initialize base, ascending, and intermediate segments as vectors of zeros
      self$trimmed_data$base <- rep(0, self$n_data - 1)
      self$trimmed_data$ascending <- rep(0, self$n_data - 1)
      self$trimmed_data$intermediate <- rep(0, self$n_data - 1)
      
      
      
    },
    
    main = function(step_coarse = list(x = 0.1, y = 0.2),
                    step_refined = list(x = 0.05, y = 0.1),
                    inflection_accuracy = 0.001) {
      self$identify_base_segments()
      print('(base 1)')
      print(self$trimmed_data$base)
      #self$check_base_profile_consistency() #TODO
      self$identify_ascending_segments()
      print('(base 2)')
      print(self$trimmed_data$base)
      self$truncate_right()  # Added method call
      print('(base 3)')
      print(self$trimmed_data$base)
      
      self$truncate_left()
      print('(base 4)')
      print(self$trimmed_data$base)
      #self$identify_intermediate_segments()
      #print('(base 4)')
      #print(self$trimmed_data$base)
      
      self$parallelogram_truncation()
      
      self$get_area_of_interest()
      
      self$find_inflection(self$roi, step_coarse$x, step_coarse$y)
      
      last <- self$result$inflection_point[1, 1]
      
      self$refine_region_of_interest()
      
      self$find_inflection(self$roi, step_refined$x, step_refined$y, fit_type = "parabolic")
      
      current <- self$result$inflection_point[1, 1]
      refine_counter <- 0
      while (abs(last - current) > inflection_accuracy) {
        #<-beyond the original paper, consider removing or making it an option.
        print(abs(last - current))
        self$refine_region_of_interest(0.5)
        self$find_inflection(self$roi, step_refined$x, step_refined$y, fit_type = "parabolic") #<-consider decreasing step size every loop
        last <- current
        current <- self$result$inflection_point[1, 1]
        refine_counter <- refine_counter + 1
        
        if (refine_counter > 1) {
          print("unstable solution")
          break
        }
      }
      print('final inflection point found:')
      print(self$result$inflection_point)
      
      self$plot_data_with_segments()  # Plotting after all processing
      
    },
    
    # Method to clean data by filtering based on the threshold on y_values
    clean_data = function() {
      x_values <- self$input_data$x
      y_values <- self$input_data$y
      above_threshold <- which(y_values > self$threshold)
      if (length(above_threshold) == 0) {
        return(list(x = x_values, y = y_values))
      }
      last_above <- max(above_threshold)
      cleaned_x <- x_values[1:last_above]
      cleaned_y <- y_values[1:last_above]
      return(list(x = cleaned_x, y = cleaned_y))
    },
    
    # Method to calculate slopes between data points
    calculate_slopes = function(data) {
      y_values <- data$y
      slopes <- numeric(length(y_values) - 1)
      for (i in 1:(length(y_values) - 1)) {
        slopes[i] <- y_values[i + 1] - y_values[i]
      }
      return(slopes)
    },
    
    # Method to identify base segments
    identify_base_segments = function() {
      slopes <- self$calculate_slopes(self$trimmed_data)
      y_values <- self$trimmed_data$y
      for (i in 1:(length(slopes))) {
        if (is.na(slopes[i]) ||
            is.na(y_values[i]) || is.na(y_values[i + 1])) {
          next  # Skip if any value is NA
        }
        if (slopes[i] <= 0 &&
            (y_values[i] <= self$threshold ||
             y_values[i + 1] <= self$threshold)) {
          self$trimmed_data$base[i] <- 1
        }
      }
      rightmost_base <- which(self$trimmed_data$base == 1)
      if (length(rightmost_base) > 0) {
        self$trimmed_data$base[1:max(rightmost_base) + 1] <- 1
      }
      print('initially identified base segment')
      print(self$trimmed_data$base)
    },
    
    # Method to check profile consistency for base segments only
    check_base_profile_consistency = function() {
      slopes <- self$calculate_slopes(self$trimmed_data)
      print(slopes)
      y_values <- self$trimmed_data$y
      for (i in which(self$trimmed_data$base == 1)) {
        if (i < length(y_values)) {
          if (is.na(y_values[i]) || is.na(y_values[i + 1])) {
            print(paste("Warning: NA value found at index", i))
            next  # Skip this iteration if NA is found
          }
          
          if (y_values[i] > self$threshold &&
              y_values[i + 1] <= self$threshold) {
            print(paste(
              "Warning: Descends across threshold at base segment index",
              i
            ))
            warning("Profile consistency at base part: no (descends across threshold)")
          }
          
          if (abs(slopes[i]) > 0.5 * self$threshold) {
            print(paste(
              "Warning: Large difference between nodes at base segment index",
              i
            ))
            warning("Profile consistency at base part: no (large difference between nodes)")
          }
        }
      }
    },
    
    # Method to identify the initial ascending segments
    identify_ascending_segments = function() {
      x_values <- self$trimmed_data$x
      y_values <- self$trimmed_data$y
      
      # Find the index where the data first crosses the threshold in an ascending manner
      start_index <- which(y_values > self$threshold &
                             c(FALSE, diff(y_values) > 0))[1]
      
      if (!is.na(start_index)) {
        # Include the segment immediately before the first crossing
        if (start_index > 1) {
          self$trimmed_data$ascending[start_index - 1] <- 1
        }
        
        # Identify the ascending segments
        for (i in start_index:(length(y_values) - 1)) {
          if (y_values[i] > self$threshold) {
            self$trimmed_data$ascending[i] <- 1
          } else {
            break
          }
        }
      }
      print('initially identified ascending segment')
      print(self$trimmed_data$ascending)
    },
    
    # Method to truncate ascending segments and update trimmed data
    truncate_right = function() {
      y_values <- self$trimmed_data$y
      x_values <- self$trimmed_data$x
      print('x_val_before')
      print(x_values)
      print('y_val_before')
      print(y_values)
      print(y_values[1])
      last_ascending_index <- max(which(self$trimmed_data$ascending == 1))
      print('last ascending index')
      print(last_ascending_index)
      print('length of y values')
      print(length(y_values))
      
      # Ensure that last_ascending_index is valid
      if (is.na(last_ascending_index) ||
          last_ascending_index == length(y_values)) {
        print('nothing to truncate!')
        return()  # Nothing to truncate
      }
      
      # Loop backwards from the end of the ascending region
      
      while (last_ascending_index > 1) {
        if (y_values[last_ascending_index + 1] > y_values[last_ascending_index]) {
          self$trimmed_data$ascending[last_ascending_index + 1] <- 1
          break  # Stop if the last segment is positively sloped
        }
        self$trimmed_data$ascending[last_ascending_index] <- 0
        last_ascending_index <- last_ascending_index - 1
      }
      
      # Truncate the trimmed data to ensure no intermediate segments exist beyond the last ascending point
      print('datarange')
      print(1:13)
      self$trimmed_data$x <- x_values[1:(last_ascending_index +
                                           1)]
      self$trimmed_data$y <- y_values[1:(last_ascending_index +
                                           1)]
      print('x_val_after')
      print(self$trimmed_data$x)
      print('y_val_after')
      print(self$trimmed_data$y)
      
      print('before truncation base')
      print(self$trimmed_data$base)
      
      print('last_index')
      print(last_ascending_index)
      self$trimmed_data$base <- self$trimmed_data$base[1:(last_ascending_index +
                                                            1)]
      print('after truncation base')
      print(self$trimmed_data$base)
      
      self$trimmed_data$ascending <- self$trimmed_data$ascending[1:(last_ascending_index +
                                                                      1)]
      self$trimmed_data$intermediate <- self$trimmed_data$intermediate[1:(last_ascending_index +
                                                                            1)]
      
      print('updated trimmed base')
      print(self$trimmed_data$base)
      
      print('updated trimmed ascending')
      print(self$trimmed_data$ascending)
      
      print('updated trimmed intermediate')
      print(self$trimmed_data$intermediate)
    },
    
    # Method to truncate ascending segments and update trimmed data
    truncate_left = function() {
      y_values <- self$trimmed_data$y
      threshold_indices <- ifelse(y_values <= self$threshold, 1, 0)
      initial_indices <- self$trimmed_data$base
      # Update the initial indices by combining with threshold results
      sequence <- initial_indices * threshold_indices
      
      start_index <- which(sequence == 1)[1]
      
      # Extract the sequence from the first 1 to the end
      self$trimmed_data$base <- self$trimmed_data$base[start_index:length(sequence)]
      print('trimmed_base_left')
      print(self$trimmed_data$base)
      self$trimmed_data$ascending <- self$trimmed_data$ascending[start_index:length(sequence)]
      print('trimmed_ascending_left')
      print(self$trimmed_data$ascending)
      
      self$trimmed_data$intermediate <- self$trimmed_data$intermediate[start_index:length(sequence)]
      
      print('trimmed_intermed_left')
      print(self$trimmed_data$intermediate)
      
      self$trimmed_data$x <- self$trimmed_data$x[start_index:length(sequence)]
      self$trimmed_data$y <- self$trimmed_data$y[start_index:length(sequence)]
      print('trimmed_left_x')
      print(self$trimmed_data$x)
      print('trimmed_left_y')
      print(self$trimmed_data$y)
      
    },
    
    # Method to identify intermediate segments
    identify_intermediate_segments = function() {
      #  print(paste("trimmed_data:", self$trimmed_data$x))
      for (i in seq_len(self$trimmed_data$x - 1)) {
        # Check if both ascending and base are not 1 at index i
        if (!is.na(self$trimmed_data$ascending[i]) &&
            !is.na(self$trimmed_data$base[i]) &&
            self$trimmed_data$ascending[i] != 1 &&
            self$trimmed_data$base[i] != 1) {
          self$trimmed_data$intermediate[i] <- 1
        }
      }
      print('identified intermediate_segments')
      print(self$trimmed_data$intermediate)
    },
    
    # Method to get data nodes from logical segment vectors
    get_nodes = function() {
      segments <- list(
        base = integer(0),
        ascending = integer(0),
        intermediate = integer(0)
      )
      for (i in seq_len(length(self$trimmed_data$x) - 1)) {
        for (key in names(segments)) {
          data_vector <- self$trimmed_data[[key]] # Changed to input_data
          
          
          # Ensure data_vector is not null and is a valid vector
          if (is.null(data_vector) ||
              !is.vector(data_vector)) {
            print(paste("Data vector for key", key, "is NULL or not a vector"))
            next
          }
          
          # Check index bounds
          if (i <= length(data_vector)) {
            if (!is.na(data_vector[i]) && data_vector[i] == 1) {
              segments[[key]] <- c(segments[[key]], i)
            }
          } else {
            
          }
        }
      }
      print('get_nodes_segments')
      print(segments)
      return(segments)
      print('get_nodes_segments')
      print(segments)
    },
    
    parallelogram_truncation = function() {
      ## TODO: ACTUALLY TRUNCATE BASED ON THIS
      x_values <- self$trimmed_data$x[which(self$trimmed_data$ascending ==
                                              1)]
      y_values <- self$trimmed_data$y[which(self$trimmed_data$ascending ==
                                              1)]
      params <- self$parallelogram_optimizer$optimize_parallelogram(x = x_values, y = y_values)
      # Plot and ensure it renders
      cat("Parallelogram Parameters:\n")
      print(params)
      self$parallelogram_params <- params
      # return(params)
    },
    
    # Method to define region of interest
    get_area_of_interest = function() {
      segments <- self$get_nodes()
      # Ensure that base and ascending segments are not empty
      if (length(segments$base) > 0 &&
          length(segments$ascending) > 0) {
        # Check for intermediate segments
        if (length(segments$intermediate) == 0) {
          # No intermediate segment
          
          # Check if the last base and first ascending segments are contiguous
          if (segments$base[length(segments$base)] == segments$ascending[1]) {
            start <- (self$trimmed_data$x[segments$base[length(segments$base)]] +
                        self$trimmed_data$x[segments$base[length(segments$base) - 1]]) / 2
            end <- 0.95 * (self$trimmed_data$x[segments$ascending[2]] -
                             self$trimmed_data$x[segments$ascending[1]]) +
              self$trimmed_data$x[segments$ascending[1]]
            
            lower <- (self$trimmed_data$y[segments$base[length(segments$base)]] +
                        self$trimmed_data$y[segments$base[length(segments$base) - 1]]) / 2
            upper <- 0.95 * (self$trimmed_data$y[segments$ascending[2]] -
                               self$trimmed_data$y[segments$ascending[1]]) +
              self$trimmed_data$y[segments$ascending[1]]
            
          } else {
            # Case 3: No nodes below the node separating base and ascending segments
            cat("i'm at case 3 !!!")
            start <- self$trimmed_data$x[segments$ascending[1]]
            end <- 0.95 * (self$trimmed_data$x[segments$ascending[2]] -
                             self$trimmed_data$x[segments$ascending[1]]) +
              self$trimmed_data$x[segments$ascending[1]]
            
            lower <- self$trimmed_data$y[segments$ascending[1]]
            upper <- 0.95 * (self$trimmed_data$y[segments$ascending[2]] -
                               self$trimmed_data$y[segments$ascending[1]]) +
              self$trimmed_data$y[segments$ascending[1]]
          }
          
        } else {
          # Case 2: Intermediate segments exist
          cat("i'm at case 2 !!!")
          # Exclude the leftmost node of the first intermediate segment
          start <- self$trimmed_data$x[segments$intermediate[2]]
          end <- 0.95 * (self$trimmed_data$x[segments$ascending[2]] -
                           self$trimmed_data$x[segments$ascending[1]]) +
            self$trimmed_data$x[segments$ascending[1]]
          
          lower <- self$trimmed_data$y[segments$intermediate[2]]
          upper <- 0.95 * (self$trimmed_data$y[segments$ascending[2]] -
                             self$trimmed_data$y[segments$ascending[1]]) +
            self$trimmed_data$y[segments$ascending[1]]
        }
        
        area <- list(x = c(start, end), y = c(lower, upper))
        
        self$roi <- area
        print('roi')
        print(self$roi)
      }
    },
    
    find_inflection = function(roi, step_x, step_y, fit_type = "linear") {
      self$result <- self$inflectionPointSeeker$seek(self$trimmed_data, roi, step_x, step_y, fit_type =
                                                       fit_type)
      print(self$roi)
      print(self$result$inflection_point[1, 1])
      print(self$result$inflection_point[1, 2])
      
    },
    
    refine_region_of_interest = function(area_reduction =
                                           0.1) {
      dx <- (self$roi$x[2] - self$roi$x[1]) * area_reduction ^ 0.5
      dy <- (self$roi$y[2] - self$roi$y[1]) * area_reduction ^
        0.5
      
      start <- self$result$inflection_point[1, 1] - dx * 0.5
      end <- self$result$inflection_point[1, 1] + dx * 0.5
      lower <- self$result$inflection_point[1, 2] - dy * 0.5
      upper <- self$result$inflection_point[1, 2] + dy * 0.5
      self$coarse_roi <- self$roi
      self$roi = list(x = c(start, end), y = c(lower, upper))
    },
    
    # Method to plot data and identified segments
    plot_data_with_segments = function() {
      x_values <- self$input_data$x
      y_values <- self$input_data$y
      trimmed_x <- self$trimmed_data$x
      trimmed_y <- self$trimmed_data$y
      
      base_seg <- self$trimmed_data$base
      base_segments <- base_seg == 1
      ascending_seg <- self$trimmed_data$ascending
      ascending_segments <- ascending_seg == 1
      
      # Create data frames for plotting
      plot_data <- data.frame(x = x_values, y = y_values)
      trimmed_plot_data <- data.frame(x = trimmed_x, y = trimmed_y)
      
      # Create shaded areas for base and ascending segments
      shaded_area_data_base <- data.frame(
        xmin = numeric(),
        xmax = numeric(),
        ymin = numeric(),
        ymax = numeric()
      )
      
      shaded_area_data_ascending <- data.frame(
        xmin = numeric(),
        xmax = numeric(),
        ymin = numeric(),
        ymax = numeric()
      )
      
      # Determine base segments area
      if (any(base_segments)) {
        first_true <- which(base_segments)[1]
        last_true <- tail(which(base_segments), 1)
        
        shaded_area_data_base <- data.frame(
          xmin = self$trimmed_data$x[first_true],
          xmax = self$trimmed_data$x[last_true],
          # Ensure right boundary
          ymin = -Inf,
          ymax = Inf
        )
        shaded_area_data_base <- shaded_area_data_base[order(shaded_area_data_base$xmin), ]
      }
      
      # Determine ascending segments area
      if (any(ascending_segments)) {
        first_true <- which(ascending_segments)[1]
        last_true <- tail(which(ascending_segments), 1)
        
        shaded_area_data_ascending <- data.frame(
          xmin = self$trimmed_data$x[first_true],
          xmax = self$trimmed_data$x[last_true],
          # Ensure right boundary
          ymin = -Inf,
          ymax = Inf
        )
        shaded_area_data_ascending <- shaded_area_data_ascending[order(shaded_area_data_ascending$xmin), ]
      }
      
      # Initialize the ggplot object
      p <- ggplot()
      
      # plot roi heat map
      if (!is.null(self$inflectionPointSeeker$grid_points)) {
        p <- p + geom_point(data = self$inflectionPointSeeker$grid_points,
                            aes(
                              x = x,
                              y = y,
                              color = log(self$inflectionPointSeeker$residuals)
                            )) +
          scale_color_gradient(low = "blue", high = "red")  # Color gradient from blue to red
        
        # overlay slope lines
        print(self$result$right_para)
        print(length(self$result$right_para))
        if (length(self$result$right_para) == 1) {
          print("fit linear right")
          p <- p + geom_abline(
            intercept = self$result$inflection_point[1, 2] - self$result$right_para * self$result$inflection_point[1, 1],
            slope = self$result$right_para,
            color = "darkgrey",
            size = 2
          )   # Right line
        }
        else{
          # curve(optimal_a * x^2 + optimal_b * x + optimal_c, from = min(x_right), to = max(x_right), add = TRUE, col = "red", lwd = 2)
          p <- p + geom_function(
            fun = function(x)
              self$result$right_para[1] * x ^ 2 + self$result$right_para[2] * x + self$result$right_para[3],
            color = "darkgrey",
            xlim = c(self$result$inflection_point[1, 1], trimmed_plot_data$x[nrow(trimmed_plot_data)]),
            size = 1.5,
            # Equivalent to lwd in base R
            n = 1000  # Number of points to evaluate for smooth curve
          )
          # coord_cartesian(xlim = c(self$result$inflection_point[1,1],trimmed_plot_data$x[nrow(trimmed_plot_data)]))
        }
        # p<-p+geom_segment(intercept = self$result$inflection_point[1,2] - self$result$left_para * self$result$inflection_point[1,1], slope = self$result$left_para, color = "darkgrey", size = 2, xlim = c(trimmed_plot_data$x[1],self$result$inflection_point[1,1])) +  # Left line
        p <- p + geom_segment(
          # m * (x - poi_x) + poi_y
          aes(
            x = trimmed_plot_data$x[1],
            y = self$result$left_para * (trimmed_plot_data$x[1] - self$result$inflection_point[1, 1]) + self$result$inflection_point[1, 2],
            xend = self$result$inflection_point[1, 1],
            yend = self$result$inflection_point[1, 2]
          ),
          color = "darkgrey",
          size = 2
        )
        p <- p + geom_point(
          aes(
            x = self$result$inflection_point[1, 1],
            y = self$result$inflection_point[1, 2]
          ),
          color = "black",
          size = 6
        )   # Plot the POI
        
      }
      
      # Raw data as black dotted line
      p <- p + geom_line(
        data = plot_data,
        aes(x = x, y = y),
        linetype = "dotted",
        color = "black",
        linewidth = 1
      ) +
        geom_point(
          data = plot_data,
          aes(x = x, y = y),
          color = "black",
          size = 2
        ) +
        
        # Trimmed data as thicker blue line
        geom_line(
          data = trimmed_plot_data,
          aes(x = x, y = y),
          color = "blue",
          linewidth = 1.5
        ) +
        geom_point(data = trimmed_plot_data,
                   aes(x = x, y = y),
                   color = "blue",
                   size = 2) +
        
        # Horizontal line for threshold
        geom_hline(
          yintercept = self$threshold,
          linetype = "dashed",
          color = "red"
        ) +
        
        # Add labels and theme
        labs(title = "Melatonin profile", x = "time of day", y = "salivary melatonin [pg/ml]") +
        theme_minimal() +
        theme(legend.position = "none")
      
      # Overlay shaded areas behind base segments and ascending segments
      if (nrow(shaded_area_data_base) > 0) {
        p <- p + geom_rect(
          data = shaded_area_data_base,
          aes(
            xmin = xmin,
            xmax = xmax,
            ymin = ymin,
            ymax = ymax
          ),
          fill = "lightblue",
          alpha = 0.3,
          color = NA
        )
      }
      if (nrow(shaded_area_data_ascending) > 0) {
        p <- p + geom_rect(
          data = shaded_area_data_ascending,
          aes(
            xmin = xmin,
            xmax = xmax,
            ymin = ymin,
            ymax = ymax
          ),
          fill = "lightpink",
          alpha = 0.3,
          color = NA
        )
      }
      
      if (!is.null(self$roi)) {
        roi_box <- data.frame(
          xmin = self$roi$x[1],
          xmax = self$roi$x[2],
          # Ensure right boundary
          ymin = self$roi$y[1],
          ymax = self$roi$y[2]
        )
        # p<- p+geom_rect(data = roi_box, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = NA, alpha = 0.3, color = "darkgreen", linewidth=2)
      }
      # Return the plot
      print(p)
    },
    
    
    # Method to plot parallelogram
    plot_parallelogram = function(existing_plot = NULL,
                                  params = NULL) {
      # Check if the optimizer is available
      if (is.null(self$parallelogram_optimizer)) {
        stop(
          "ParallelogramOptimizer instance is not available. Run parallelogram_truncation() first."
        )
      }
      # Plot the parallelogram using the stored optimizer
      # self$parallelogram_optimizer$plot_parallelogram(existing_plot, params)
    }
  )
)

# Example usage:
# all_data <- list(
#   x = c(19, 19.5, 20, 20.5, 21, 21.5, 22, 22.5, 23, 23.5, 24, 24.5, 25, 25.5, 26),
#   y = c(5, 0.1, 2.0, 1.7, 0.8, 1.1, 1, 5, 9, 15, 36, 36.5, 38, 22, 1)
# )
all_data <- list(
  # x = c(19, 19.5, 20, 20.5, 21, 21.5, 22, 22.5, 23, 23.5, 24, 24.5, 25),
  # y = c(5, 0.1, 2.0, 1.7, 0.8, 1.1, 1, 21, 23, 25, 36, 36.5, 38)
  #x = c(19, 19.5, 20, 20.5, 21, 21.5, 22, 22.5, 23, 23.5, 24, 24.5, 25), #
  #y = c(5, 0.1, 2.0, 1.7, 0.8, 1.1, 1, 1.1, 3, 5, 36, 36.5, 38)  # this broke the base assignment
  x = c(18, 18.5, 19, 19.5, 20, 20.5, 21, 21.5, 22, 22.5, 23.5, 24, 25, 26, 27),
  y = c(6, 5, 2, 0.1, 2.0, 1.7, 0.8, 1.1, 1, 5, 10, 20, 38, 15, 1)
  
)

threshold <- 2.3

# Create a DataProcessor object
processor <- DLMO$new(all_data, threshold)
processor$main()

# Print results
# print("input data:")
# print(processor$input_data)
# print("Base segments:")
# print(processor$input_data$base)
# print("Ascending segments:")
# print(processor$input_data$ascending)
# print("Intermediate segments:")
# print(processor$input_data$intermediate)

#cleanedupdatax = processor$trimmed_data$x
#cleanedupdatay = processor$trimmed_data$y

#print(paste("trim x:", cleanedupdatax))
#print(paste("trim y:", cleanedupdatay))
# Plot the raw data with base and ascending segments highlighted
base_plot <- processor$plot_data_with_segments()

# Perform parallelogram truncation on the truncated ascending segments and plot the parallelogram

# overlay parallelogram plot

processor$plot_parallelogram(existing_plot = base_plot,
                             params = processor$parallelogram_params)

plot


x = c(18, 18.5, 19, 19.5, 20, 20.5, 21, 21.5, 22, 22.5, 23.5, 24, 25)
y = c(6, 5, 2, 0.1, 2.0, 1.7, 0.8, 1.1, 1, 5, 10, 20, 38)

# Create a data frame
df <- data.frame(x = x, y = y)

# Plot using ggplot2
library(ggplot2)
ggplot(df, aes(x = x, y = y)) +
  geom_line(color = "blue") +
  geom_point() +
  labs(title = "Y vs Indices", x = "Index", y = "Y Value")
