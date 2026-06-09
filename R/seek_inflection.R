#' Generate a Grid of Points Within a Region of Interest
#'
#' This function generates a uniform grid of points within a specified region of interest (ROI).
#' The grid is created based on user-defined step sizes for the x and y dimensions.
#'
#' @param roi A list containing:
#'   - `x`: A numeric vector of length 2 specifying the x-axis range (`xmin`, `xmax`).
#'   - `y`: A numeric vector of length 2 specifying the y-axis range (`ymin`, `ymax`).
#' @param step_x Numeric. The step size for x-coordinates (grid resolution along the x-axis).
#' @param step_y Numeric. The step size for y-coordinates (grid resolution along the y-axis).
#'
#' @return A data frame containing all grid points with columns:
#'   - `x`: The x-coordinates of the grid points.
#'   - `y`: The y-coordinates of the grid points.
#'
#' @seealso \code{\link{seek_inflection}}, which utilizes this function to generate search grids.
#'
#' @examples
#' roi <- list(x = c(1, 5), y = c(2, 6))
#' grid <- make_grid(roi, step_x = 0.5, step_y = 0.5)
#' print(grid)
#'
#' @export
make_grid <- function(roi, step_x, step_y) {
  # Extract the bounding limits of the region of interest (ROI)
  xmin <- roi$x[1]  # Minimum x-value
  xmax <- roi$x[2]  # Maximum x-value
  ymin <- roi$y[1]  # Minimum y-value
  ymax <- roi$y[2]  # Maximum y-value

  if (!all(is.finite(c(xmin, xmax, ymin, ymax)))) {
    stop(
      "Cannot create the DLMO search grid because the region of interest contains non-finite bounds. ",
      "This usually means the base, intermediate, or ascending segments were not valid for the selected ",
      "`threshold` and `interval_limit`. Inspect the segmented profile and consider increasing ",
      "`interval_limit` to include a later sustained rise, or adjusting `threshold` so the selected rise ",
      "better matches the intended DLMO event.",
      call. = FALSE
    )
  }

  if (xmin >= xmax || ymin >= ymax) {
    stop(
      "Cannot create the DLMO search grid because the region of interest has invalid bounds. ",
      "The start must be before the end, and the lower melatonin bound must be below the threshold. ",
      "Inspect the segmented profile and consider adjusting `threshold` and/or `interval_limit`.",
      call. = FALSE
    )
  }

  # Generate sequences of x and y values based on step sizes
  x_seq <- seq(from = xmin, to = xmax, by = step_x)  # X-coordinates
  y_seq <- seq(from = ymin, to = ymax, by = step_y)  # Y-coordinates

  # Create a data frame containing all combinations of x and y values
  grid_points <- expand.grid(x = x_seq, y = y_seq)

  # Return the generated grid of points
  return(grid_points)
}


#' Nonlinear Constraints for Parabolic Fitting
#'
#' This function enforces constraints on a parabolic fit for ascending melatonin concentration data.
#' The constraints ensure that the fitted parabola passes through a specified point of interest (POI)
#' and maintains a positive slope in the ascending segment.
#'
#' @param params A numeric vector of length 3 representing the parabola parameters:
#'   - `a`: Quadratic coefficient.
#'   - `b`: Linear coefficient.
#'   - `c`: Constant term.
#' @param poi_x Numeric. The x-coordinate of the point of interest (POI).
#' @param poi_y Numeric. The y-coordinate of the POI.
#' @param x Numeric vector of x-values representing the ascending segment.
#'
#' @return A list containing:
#'   - `pass_poi`: The squared deviation from the constraint that the parabola must pass through the POI.
#'   - `pos_grad`: The squared deviation from the constraint ensuring a positive slope at the last ascending point.
#'
#' @details
#' The function imposes two constraints:
#' 1. **Passing through POI**: Ensures the fitted parabola satisfies \eqn{y = ax^2 + bx + c} at `poi_x, poi_y`.
#' 2. **Positive Slope**: Ensures \eqn{\frac{dy}{dx} = 2ax + b} remains non-negative at the last x-value.
#'
#' If the minimum slope at the last ascending point is negative, it is forced to zero as a penalty.
#'
#' @seealso \code{\link{objective_function}}, which utilizes these constraints during optimization.
#'
#' @examples
#' params <- c(0.1, 2, 3)  # Example parabola: y = 0.1x^2 + 2x + 3
#' poi_x <- 5
#' poi_y <- 15
#' x_values <- seq(4, 6, by = 0.5)
#' constraints <- nl_constraints(params, poi_x, poi_y, x_values)
#' print(constraints)
#'
#' @export
nl_constraints <- function(params, poi_x, poi_y, x) {
  # Extract quadratic, linear, and constant terms from parameters
  a <- params[1]  # Quadratic coefficient
  b <- params[2]  # Linear coefficient
  c <- params[3]  # Constant term

  # Constraint 1: Ensure the parabola passes through the Point of Interest (POI)
  poi_diff <- a * poi_x^2 + b * poi_x + c - poi_y  # Difference from expected y-value

  # Constraint 2: Ensure the slope (dy/dx) remains positive at the last ascending x-value
  slope_values <- 2 * a * x + b  # Compute derivative of y = ax^2 + bx + c
  slope_min <- tail(slope_values, n = 1)  # Evaluate slope at the last x-value

  # If the minimum slope is negative, force it to zero as a penalty
  if (slope_min >= 0) {
    slope_min <- 0
  }

  # Return squared penalties for constraint violations
  return(list(
    pass_poi = poi_diff^2,  # Squared deviation from POI constraint
    pos_grad = slope_min^2   # Squared deviation for positive slope constraint
  ))
}


#' Enforce Slope Constraints at the Point of Interest (POI)
#'
#' This function ensures that the slope at the Point of Interest (POI) falls within the specified bounds.
#' It applies different calculations depending on whether the fit is linear or parabolic.
#'
#' @param params A numeric vector containing the model parameters:
#'   - For a linear fit: `params[1]` is the slope.
#'   - For a parabolic fit: `params[1]` (quadratic coefficient) and `params[2]` (linear coefficient).
#' @param slope_bounds A numeric vector of length 2 specifying the lower and upper slope limits.
#' @param poi_x Numeric. The x-coordinate of the Point of Interest.
#' @param fit_type Character. Specifies the type of fit used. Can be `"linear"` or `"parabolic"`.
#'
#' @return A numeric value:
#'   - `0` if the computed slope is within bounds.
#'   - A squared penalty value if the computed slope is outside the allowable range.
#'
#' @details
#' - If the fit type is `"linear"`, the slope is taken directly from `params[1]`.
#' - If the fit type is `"parabolic"`, the slope at the POI is computed as \eqn{2a \cdot poi_x + b}.
#' - The function checks if the computed slope falls within the specified `slope_bounds`.
#' - If the slope is out of bounds, it applies a squared penalty to encourage corrections during optimization.
#'
#' @seealso \code{\link{objective_function}} which utilizes this constraint during optimization.
#'
#' @examples
#' params_linear <- c(0.5)  # Example slope for a linear fit
#' slope_bounds <- c(0.2, 1.0)
#' poi_x <- 5
#' poi_constraints(params_linear, slope_bounds, poi_x, fit_type = "linear")
#'
#' params_parabolic <- c(0.1, 0.5)  # Example quadratic and linear coefficients
#' poi_constraints(params_parabolic, slope_bounds, poi_x, fit_type = "parabolic")
#'
#' @export
poi_constraints <- function(params, slope_bounds, poi_x, fit_type = "linear") {
  # Determine the slope value based on fit type
  if (fit_type == "linear") {
    slope_value <- params[1]  # Directly extract slope for a linear fit
  } else {  # Parabolic fit
    a <- params[1]  # Quadratic coefficient
    b <- params[2]  # Linear coefficient
    slope_value <- 2 * a * poi_x + b  # Compute slope at POI for a parabolic function
  }

  # Constraint: Ensure slope at POI falls within the allowed range
  if (dplyr::between(slope_value, slope_bounds[1], slope_bounds[2])) {
    return(0)  # No penalty if within bounds
  } else {
    # Compute squared penalty for out-of-bounds slopes
    return(min((slope_bounds[1] - slope_value)^2, (slope_bounds[2] - slope_value)^2))
  }
}


#' Constraint Function for Base Fit Edge Alignment
#'
#' This function ensures that the base fit remains within defined vertical bounds
#' by applying penalties if the left or right edge of the base segment exceeds limits.
#'
#' @param y Numeric vector representing the fitted y-values along the base segment.
#' @param min_y Numeric. The minimum y-value in the dataset (used as a lower bound).
#' @param poi_y Numeric. The y-value at the Point of Interest (POI).
#' @param threshold Numeric. The upper y-bound for the base fit. Default is `2.3`.
#'
#' @return A numeric penalty value:
#'   - `0` if both edges are within the defined range.
#'   - A squared penalty if the left or right edges exceed their respective bounds.
#'
#' @details
#' - The function applies constraints separately for the left and right edges of the base segment.
#' - The left constraint ensures that the starting y-value remains within `[0, threshold]`.
#' - The right constraint ensures that the POI y-value is between `[min_y, threshold]`.
#' - If either constraint is violated, a squared penalty is applied to encourage corrections during optimization.
#'
#' @seealso \code{\link{objective_function}} where this function is used in the optimization process.
#'
#' @examples
#' y_values <- c(0.5, 1.2, 2.0)  # Example y-values for the base segment
#' min_y <- 0.1
#' poi_y <- 1.8
#' threshold <- 2.3
#' base_constraint(y_values, min_y, poi_y, threshold)
#'
#' @export
base_constraint <- function(y, min_y, poi_y, threshold = 2.3) {
  # Ensure the left edge of the base segment is within the allowed range [0, threshold]
  if (dplyr::between(y[1], 0, threshold)) {
    left_constr <- 0  # No penalty if within bounds
  } else {
    # Apply squared penalty if outside bounds
    left_constr <- min((y[1] - threshold)^2, y[1]^2)
  }

  # Ensure the right edge at POI is between [min_y, threshold]
  if (dplyr::between(poi_y, min_y, threshold)) {
    right_constr <- 0  # No penalty if within bounds
  } else {
    # Apply squared penalty if outside bounds
    right_constr <- min((poi_y - threshold)^2, (poi_y - min_y)^2)
  }

  # Return total penalty (sum of left and right constraint violations)
  return(left_constr + right_constr)
}



#' Objective Function for Line or Parabola Fitting
#'
#' This function calculates the cost for fitting either a linear or parabolic model
#' to a segment of the melatonin profile. The objective minimizes the residual error
#' while ensuring the fit adheres to predefined constraints.
#'
#' @param params Numeric vector of parameters to optimize:
#'   - If `fit_type == "linear"`, `params[1]` is the slope (m).
#'   - If `fit_type == "parabolic"`, `params[1:2]` represent coefficients (a, b), and c is computed.
#' @param x Numeric vector of x-coordinates (time points in decimal format).
#' @param y Numeric vector of y-coordinates (melatonin concentrations).
#' @param poi List containing:
#'   - `x`: x-coordinate of the Point of Interest (POI).
#'   - `y`: y-coordinate of the POI.
#' @param fit_type Character. Either `"linear"` or `"parabolic"` to specify the fitting model.
#' @param slope_bounds Numeric vector of length 2 defining the allowed slope range.
#' @param base_id Numeric vector indicating which points belong to the base segment (1 = base, 0 = ascending).
#' @param region Character. Either `"base"` or `"ascending"` to specify the segment being fitted.
#' @param weight_base Logical. If `TRUE`, assigns a lower weight to base points in the residual computation.
#' @param threshold Numeric. Upper bound constraint for base fitting. Default is `threshold`.
#'
#' @return Numeric value representing the total cost (sum of squared residuals and constraint penalties).
#'
#' @details
#' - The function supports both linear and parabolic fits.
#' - The `base_constraint` function penalizes fits that exceed the threshold in the base region.
#' - The `nl_constraints` function ensures that parabolic fits maintain an increasing slope.
#' - The `poi_constraints` function ensures that the fit slope at the POI stays within the given bounds.
#' - Weighted residuals are applied to emphasize ascending points in parabolic fitting.
#'
#' @seealso
#' \code{\link{base_constraint}}, \code{\link{nl_constraints}}, \code{\link{poi_constraints}}
#'
#' @examples
#' x <- c(0, 1, 2, 3, 4)
#' y <- c(0.5, 1.0, 2.0, 3.5, 5.0)
#' poi <- list(x = 2, y = 2.0)
#' params <- c(1.5)  # Linear slope example
#' slope_bounds <- c(-2, 2)
#' objective_function(params, x, y, poi, "linear", slope_bounds, base_id = c(1,1,0,0,0), "ascending", weight_base = TRUE)
#'
#' @export
objective_function <- function(params, x, y, poi, fit_type, slope_bounds, base_id, region, weight_base, threshold = threshold) {
  poi_x <- poi$x
  poi_y <- poi$y

  # Linear fitting case
  if (fit_type == "linear" & region == "ascending") {
    m <- params[1]  # Slope of the linear fit
    y_pred <- m * (x - poi_x) + poi_y
    constr_cost <- 0  # No additional constraints for ascending linear fits

  } else if (fit_type == "linear") {  # Base segment linear fit
    m <- params[1]
    y_pred <- m * (x - poi_x) + poi_y
    # Apply constraint penalty if base segment exceeds allowed bounds
    constr_cost <- base_constraint(y_pred, min(y), poi_y, threshold = threshold) * 100

  } else {  # Parabolic fitting case
    a <- params[1]
    b <- params[2]
    c <- poi$y - a * poi$x^2 - b * poi$x  # Compute c explicitly
    params <- c(a, b, c)

    # Apply nonlinear constraints to ensure the parabola meets required conditions
    constraints <- nl_constraints(params, poi_x, poi_y, x)
    constr_cost <- 100 * constraints$pos_grad  # Penalize non-monotonicity

    # Compute predicted y-values using the quadratic equation
    y_pred <- a * x^2 + b * x + c
  }

  # Apply constraint to enforce slope bounds at POI
  poi_cost <- poi_constraints(params, slope_bounds, poi_x, fit_type)
  constr_cost <- constr_cost + poi_cost * 100  # Scale constraint penalty

  # Compute residuals (difference between actual and predicted values)
  residuals <- y - y_pred

  # Apply weighting: Give ascending points higher importance in the parabolic fit
  if (weight_base) {
    residuals <- residuals * (1 - base_id) + 0.5 * residuals * base_id
  }

  # Compute final loss: Mean squared error + constraint penalties
  l2_cost <- mean(residuals^2)

  return(l2_cost + constr_cost)  # Return total cost to be minimized
}


#' Fit a Melatonin Profile to Base or Ascending Segment
#'
#' This function fits a linear or parabolic model to either the **base** or **ascending**
#' segment of a melatonin profile relative to a specified **point of interest (POI)**.
#' It ensures that the fit adheres to predefined constraints.
#'
#' @param x Numeric vector of x-coordinates (decimal time).
#' @param y Numeric vector of y-coordinates (melatonin concentration).
#' @param poi List containing:
#'   - `x`: The x-coordinate of the Point of Interest (POI).
#'   - `y`: The y-coordinate of the POI.
#' @param fit_type Character. `"linear"` (default) or `"parabolic"` for the type of fit.
#' @param region Character. `"base"` (default) or `"ascending"` to specify which part of the profile to fit.
#' @param base_tangent Numeric. Estimated slope of the base segment (used if `fit_type = "parabolic"`).
#' @param ascending_linear_tangent Numeric. Estimated slope of the ascending segment (used if `fit_type = "parabolic"`).
#' @param base_id Numeric vector. Indicator variable where `1` represents base points and `0` represents ascending points.
#' @param weight_base Logical. If `TRUE`, applies lower weights to base points in residual computation.
#' @param threshold Numeric. Upper boundary for the base fit. Default is `threshold`.
#'
#' @return A list containing:
#'   - `residual`: The mean squared residual error.
#'   - `params`: The fitted model parameters.
#'
#' @details
#' - **Base segment fitting (`region = "base"`)**: Uses a **linear** fit constrained within slope bounds.
#' - **Ascending segment fitting (`region = "ascending"`)**:
#'   - If `fit_type = "linear"`, a **linear fit** is applied with no slope constraints.
#'   - If `fit_type = "parabolic"`, a **quadratic fit** is optimized while enforcing **monotonicity**.
#' - The function **penalizes** fits that do not meet slope constraints.
#' - **Weighted residuals** prioritize **ascending points** when computing the fit error.
#'
#' @seealso
#' \code{\link{fit_linear}}, \code{\link{objective_function}}
#'
#' @examples
#' x <- c(0, 1, 2, 3, 4)
#' y <- c(0.5, 1.0, 2.0, 3.5, 5.0)
#' poi <- list(x = 2, y = 2.0)
#' fit_profile(x, y, poi, fit_type = "linear", region = "ascending", weight_base = TRUE)
#'
#' @export
fit_profile <- function(x, y, poi, fit_type = "linear", region = "base", base_tangent = NULL, ascending_linear_tangent = NULL, base_id = NULL, weight_base = TRUE, threshold = threshold) {

  if (region == "base") {
    # Define slope bounds for base fit (ensures gentle slope)
    slope_bounds <- c(-0.2, 0.2)
    edge_bounds <- list(left = c(0, threshold), right = c(min(y), threshold))

    # Perform linear fit for base segment
    result <- fit_linear(x, y, poi, base_id, slope_bounds, edge_bounds, weight_base = weight_base)
    params <- result$params
    y_pred <- params[1] * (x - poi$x) + poi$y  # Linear equation: y = m(x - x_POI) + y_POI

  } else if (region == "ascending" & fit_type == "linear") {
    # Define slope bounds for ascending linear fit (no restrictions)
    slope_bounds <- c(-Inf, Inf)

    # Perform linear fit for ascending segment
    result <- fit_linear(x, y, poi, base_id, slope_bounds, weight_base = weight_base)
    params <- result$params
    y_pred <- params[1] * (x - poi$x) + poi$y

  } else {  # Parabolic fit for the ascending segment
    # Set slope constraints to ensure a **monotonically increasing parabola**
    slope_lowerbound <- max(c(0, base_tangent, 0.5 * ascending_linear_tangent))
    slope_upperbound <- Inf
    slope_bounds <- c(slope_lowerbound, slope_upperbound)

    # Initialize parameters for parabolic optimization (start with linear estimate)
    initial_params <- c(a = 0, b = ascending_linear_tangent)

    # Optimize the parabolic fit using L-BFGS-B method
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
      weight_base = weight_base,
      threshold = threshold,
      method = "L-BFGS-B"
    )

    # Extract optimized parameters for quadratic equation: y = ax² + bx + c
    a <- optim_result$par[1]
    b <- optim_result$par[2]
    c <- poi$y - a * poi$x^2 - b * poi$x  # Solve for c using POI

    # Compute predicted y-values
    y_pred <- a * x^2 + b * x + c
    params <- list(a = a, b = b, c = c)
  }

  # Compute residuals (difference between actual and predicted values)
  delta <- y_pred - y

  # Apply weighting: prioritize ascending points in residual computation
  if (weight_base) {
    residual <- delta * (1 - base_id) + 0.5 * delta * base_id
  } else {
    residual <- delta
  }

  # Return mean squared residual error and fitted parameters
  return(list(residual = mean(residual^2), params = params))
}


#' Perform a Weighted Linear Fit
#'
#' This function fits a **linear model** to a given segment of a melatonin profile.
#' It estimates the slope of the best-fit line while enforcing **slope constraints**
#' and **edge boundary conditions** if specified.
#'
#' @param x Numeric vector of x-coordinates (decimal time).
#' @param y Numeric vector of y-coordinates (melatonin concentration).
#' @param poi List containing:
#'   - `x`: The x-coordinate of the Point of Interest (POI).
#'   - `y`: The y-coordinate of the POI.
#' @param base_id Numeric vector. Indicator variable where `1` represents base points and `0` represents ascending points.
#' @param slope_bounds Numeric vector of length 2 specifying lower and upper bounds for the slope. Default is `NULL` (no constraints).
#' @param edge_bounds List containing:
#'   - `left`: Numeric vector `[min, max]` range for left boundary (base segment).
#'   - `right`: Numeric vector `[min, max]` range for right boundary (ascending segment).
#'   If `NULL`, no edge constraints are applied.
#' @param weight_base Logical. If `TRUE`, applies lower weights to base points in residual computation.
#'
#' @return A list containing:
#'   - `params`: A numeric vector with the estimated slope.
#'
#' @details
#' - The function **applies weights** to prioritize ascending points over base points.
#' - **Slope constraints** ensure that the estimated slope remains within predefined bounds.
#' - If **edge constraints** are provided, the function ensures that the fitted line **does not exceed** specified y-limits at the segment edges.
#'
#' @seealso
#' \code{\link{fit_profile}}, \code{\link{objective_function}}
#'
#' @examples
#' x <- c(0, 1, 2, 3, 4)
#' y <- c(0.5, 1.0, 2.0, 3.5, 5.0)
#' poi <- list(x = 2, y = 2.0)
#' base_id <- c(1, 1, 0, 0, 0)
#' fit_linear(x, y, poi, base_id = base_id, slope_bounds = c(-0.2, 0.2), weight_base = TRUE)
#'
#' @export
fit_linear <- function(x, y, poi, base_id = NULL, slope_bounds = NULL, edge_bounds = NULL, weight_base = TRUE) {

  # Compute differences relative to POI
  x_diff <- x - poi$x
  y_diff <- y - poi$y

  # Apply **weighting**: prioritize ascending points (lower weight for base points)
  if (weight_base) {
    if (is.null(base_id)) {
      w <- 1  # Uniform weighting (fallback if base_id not provided)
    } else {
      w <- (1 - base_id) + 0.5 * base_id
    }
  } else {
    w <- 1  # Uniform weighting
  }

  # Compute weighted **least squares slope**: m = Σ(w * Δx * Δy) / Σ(w * Δx²)
  num <- sum(w * x_diff * y_diff)   # Weighted numerator (covariance)
  denom <- sum(w * x_diff^2)        # Weighted denominator (variance)

  # Guard against degenerate denominators (e.g., all x equal to poi$x)
  if (!is.finite(denom) || denom == 0) {
    slope <- 0
  } else {
    slope <- num / denom
    if (!is.finite(slope)) slope <- 0
  }

  # **Enforce slope constraints**: ensure slope remains within [min, max] bounds
  if (!is.null(slope_bounds) && length(slope_bounds) == 2 && all(is.finite(slope_bounds))) {
    if (slope < min(slope_bounds)) {
      slope <- min(slope_bounds)
    } else if (slope > max(slope_bounds)) {
      slope <- max(slope_bounds)
    }
  }

  # **Apply edge constraints** (if specified)
  if (!is.null(edge_bounds)) {

    # Compute predicted left edge (x[1] is the first data point in the segment)
    left_edge <- slope * (x[1] - poi$x) + poi$y
    right_edge <- poi$y  # Right edge is fixed at the POI

    # **Enforce left edge constraints**
    if (!dplyr::between(left_edge, edge_bounds$left[1], edge_bounds$left[2])) {
      if (left_edge < edge_bounds$left[1]) {
        left_edge <- edge_bounds$left[1]  # Clip to lower bound
      } else {
        left_edge <- edge_bounds$left[2]  # Clip to upper bound
      }
      if ((x[1] - poi$x) != 0) {
        slope <- (left_edge - right_edge) / (x[1] - poi$x)  # Recompute slope
      } else {
        slope <- 0
      }
    }

    # **Enforce right edge constraints**
    if (!dplyr::between(right_edge, edge_bounds$right[1], edge_bounds$right[2])) {
      if (right_edge < edge_bounds$right[1]) {
        right_edge <- edge_bounds$right[1]  # Clip to lower bound
      } else {
        right_edge <- edge_bounds$right[2]  # Clip to upper bound
      }
      if ((x[1] - poi$x) != 0) {
        slope <- (left_edge - right_edge) / (x[1] - poi$x)  # Recompute slope
      } else {
        slope <- 0
      }
    }
  }

  # Return the estimated slope
  return(list(params = c(unname(slope))))
}

#' Fit Two Splines Around a Point of Interest (POI)
#'
#' This function fits a **piecewise linear or parabolic model** to melatonin profile data.
#' It separately models the **base segment** and the **ascending segment**, ensuring a smooth transition at the **point of interest (POI)**.
#'
#' @param data A tibble containing melatonin concentration data with `datetime` values.
#' @param poi A list containing:
#'   - `x`: The x-coordinate (decimal time) of the POI.
#'   - `y`: The y-coordinate (melatonin concentration) of the POI.
#' @param fit_type Character. Type of fit for the ascending segment:
#'   - `"linear"` (default): Uses a linear model for both base and ascending segments.
#'   - `"parabolic"`: Uses a **parabolic fit** for the ascending segment.
#' @param threshold Numeric. The threshold for constraining the base segment.
#'
#' @return A list containing:
#'   - `residual`: The sum of squared residuals from both fits.
#'   - `base_params`: The estimated parameters for the base fit.
#'   - `ascending_params`: The estimated parameters for the ascending fit.
#'
#' @details
#' - The **base segment** is always fit using a **linear model**.
#' - The **ascending segment** is fit using a **linear** or **parabolic** model depending on `fit_type`.
#' - If the **ascending segment** has at least **3 data points**, it first fits a linear model and **then refines it** with a parabolic fit.
#' - This function is used in **inflection point detection**.
#'
#' @seealso
#' \code{\link{fit_profile}}, \code{\link{objective_function}}, \code{\link{seek_inflection}}
#'
#' @examples
#' data <- tibble::tibble(
#'   datetime = as.POSIXct(c("2024-06-01 22:00", "2024-06-01 22:15", "2024-06-01 22:30")),
#'   melatonin = c(0.5, 0.8, 1.2),
#'   base = c(1, 1, 0)  # Base segment indicators
#' )
#' poi <- list(x = posixct_to_decimal(data$datetime[2], data$datetime), y = data$melatonin[2])
#' fit(data, poi, fit_type = "linear", threshold = 2.3)
#'
#' @export
fit <- function(data, poi, fit_type = "linear", threshold = threshold) {

  # Convert datetime to decimal time (reference to dataset)
  x <- posixct_to_decimal(data$datetime, data$datetime[3])
  y <- data$melatonin
  base_id <- data$base  # Indicator for base segment

  # Extract POI coordinates
  poi_x <- poi$x
  poi_y <- poi$y

  # Disable weighting for base points
  weight_base <- FALSE

  ### Step 1: Fit Base Segment (Left of POI)
  left_indcs <- which(x <= poi_x)  # Identify points before or at the POI
  result_base <- fit_profile(
    x = x[left_indcs], y = y[left_indcs], poi = poi,
    fit_type = "linear", region = "base",
    base_id = base_id[left_indcs], weight_base = weight_base,
    threshold = threshold
  )

  ### Step 2: Fit Ascending Segment (Right of POI)
  right_indcs <- which(x > poi_x)  # Identify points after the POI
  result_ascending <- fit_profile(
    x = x[right_indcs], y = y[right_indcs], poi = poi,
    fit_type = "linear", region = "ascending",
    base_id = base_id[right_indcs], weight_base = weight_base,
    threshold = threshold
  )

  # Optional parabolic refinement ONLY when requested (refinement phase)
  if (fit_type == "parabolic" && length(right_indcs) >= 2) {
    slope_initial_ascending <- result_ascending$params[1]  # Initial slope from linear fit
    result_ascending <- fit_profile(
      x = x[right_indcs], y = y[right_indcs], poi = poi,
      fit_type = "parabolic",
      base_tangent = -result_base$params[1],  # Use base segment slope for continuity
      ascending_linear_tangent = slope_initial_ascending,
      region = "ascending", base_id = base_id[right_indcs],
      weight_base = weight_base, threshold = threshold
    )
  }

  # Compute **total residuals** as the sum of base and ascending residuals
  total_residuals <- result_base$residual + result_ascending$residual

  # Return results
  return(list(
    residual = total_residuals,
    base_params = result_base$params,
    ascending_params = result_ascending$params
  ))
}


#' Seek the Inflection Point in a Melatonin Profile
#'
#' This function identifies the **inflection point** in a melatonin profile by performing a
#' **grid search** over a region of interest (ROI). It first conducts a **coarse search** over
#' a large grid and, if `fine_flag = TRUE`, refines the search using a **fine-grid search**.
#'
#' @param data A tibble containing melatonin profile data.
#' @param threshold Numeric. The threshold value for the fit constraint.
#' @param roi A list defining the **region of interest** (ROI) for the search:
#'   - `x`: A numeric vector of two values representing the min and max x-coordinates.
#'   - `y`: A numeric vector of two values representing the min and max y-coordinates.
#' @param step_x Numeric. The step size for x-coordinates in the **coarse grid** (default = 0.1).
#' @param step_y Numeric. The step size for y-coordinates in the **coarse grid** (default = 0.2).
#' @param step_size_small Numeric. The step size for **fine grid refinement** (default = 0.01).
#' @param fit_type Character. Type of fit used in `fit()`:
#'   - `"linear"` (default): Uses a linear model.
#'   - `"parabolic"`: Uses a parabolic model for the ascending segment.
#' @param fine_flag Logical. If `TRUE`, performs an additional **fine-grid search** (default = `fine_flag`).
#'
#' @return A list containing:
#'   - `inflection_point`: The estimated inflection point (`x`, `y`).
#'   - `base_params`: Parameters of the **base segment** fit.
#'   - `ascending_params`: Parameters of the **ascending segment** fit.
#'   - `grid_big`: The coarse grid used in the first search.
#'   - `res_big`: The residuals corresponding to `grid_big`.
#'   - `grid_small`: The fine grid used in the second search (if `fine_flag = TRUE`).
#'   - `res_small`: The residuals corresponding to `grid_small` (if `fine_flag = TRUE`).
#'
#' @details
#' - The **coarse grid search** is performed over the ROI.
#' - If `fine_flag = TRUE`, the **best 10% of coarse-grid points** are used to refine the search.
#' - The inflection point is selected as the point that **minimizes the residual error** from `fit()`.
#'
#' @seealso
#' \code{\link{make_grid}}, \code{\link{reduce_grid}}, \code{\link{fit}}
#'
#' @export
seek_inflection <- function(data, threshold = threshold, roi, step_x = 0.1, step_y = 0.2, step_size_small = 0.01, fit_type = "linear", fine_flag = fine_flag) {

  # Generate a **coarse grid** of points within the Region of Interest (ROI)
  grid_points <- make_grid(roi, step_x, step_y)

  # Initialize variables to track the best fit
  best_residual_coarse <- Inf
  best_point_coarse <- NULL
  best_params_base_coarse <- NULL
  best_params_ascending_coarse <- NULL

  best_residual_fine <- Inf
  best_point_fine <- NULL
  best_params_base_fine <- NULL
  best_params_ascending_fine <- NULL

  # Store coarse grid results
  res_big <- numeric(nrow(grid_points))
  grid_big <- grid_points

  ### Step 1: Perform **Coarse Grid Search**
  for (i in seq_len(nrow(grid_points))) {
    poi <- grid_points[i, ]  # Evaluate each grid point
    result <- fit(data, poi, fit_type = "linear", threshold = threshold)
    res_big[i] <- result$residual

    # If the new residual is smaller, update the best fit
    if (is.finite(result$residual) && result$residual < best_residual_coarse) {
      best_residual_coarse <- result$residual
      best_point_coarse <- poi
      best_params_base_coarse <- result$base_params
      best_params_ascending_coarse <- result$ascending_params
    }
  }

  if (is.null(best_point_coarse)) {
    stop("No valid coarse DLMO fit was found for this profile.")
  }

  ### Step 2: Define Function to **Refine Grid Based on Best 10%**
  reduce_grid <- function(res_big, grid_big, step_x, step_y, step_size_small, threshold = threshold) {

    # Select **top 10% of grid points** with lowest residuals
    best_10_per <- order(res_big)[1:as.integer(0.1 * nrow(grid_big))]
    best_points <- grid_big[best_10_per, ]

    # Define new ROI bounds based on selected points
    min_y <- max(min(best_points$y) - step_y, min(grid_big$y))
    max_y <- min(max(best_points$y) + step_y, max(grid_big$y))
    min_x <- max(min(best_points$x) - step_x, min(grid_big$x))
    max_x <- min(max(best_points$x) + step_x, max(grid_big$x))

    # Generate **fine grid** over new ROI
    roi <- list(x = c(min_x, max_x), y = c(min_y, max_y))
    return(make_grid(roi, step_size_small, step_size_small))
  }

  # Apply **grid refinement**
  grid_points <- reduce_grid(res_big, grid_big, step_x, step_y, step_size_small, threshold)
  grid_small <- grid_points
  res_small <- numeric(nrow(grid_points))

  ### Step 3: Perform **Fine Grid Search** (if `fine_flag = TRUE`)
  if (fine_flag) {
    for (i in seq_len(nrow(grid_points))) {
      poi <- grid_points[i, ]  # Evaluate each grid point
      result <- fit(data, poi, fit_type = "parabolic", threshold = threshold)
      res_small[i] <- result$residual

      # If the new residual is smaller, update the best fit
      if (is.finite(result$residual) && result$residual < best_residual_fine) {
        best_residual_fine <- result$residual
        best_point_fine <- poi
        best_params_base_fine <- result$base_params
        best_params_ascending_fine <- result$ascending_params
      }
    }

    if (is.null(best_point_fine)) {
      warning("No valid fine DLMO fit was found; returning coarse fit only.")
    }
  }

  ### Step 4: Return Best Fit Results
  return(list(
    inflection_point_coarse = best_point_coarse,
    inflection_point_fine = best_point_fine,
    base_params_coarse = best_params_base_coarse,
    ascending_params_coarse = best_params_ascending_coarse,
    base_params_fine = best_params_base_fine,
    ascending_params_fine = best_params_ascending_fine,
    grid_big = grid_big,
    res_big = res_big,
    grid_small = grid_small,
    res_small = res_small,
    datetime_ref = attr(roi, "datetime_ref")
  ))
}

#' Identify the Inflection Point in a Melatonin Profile
#'
#' This function determines the **inflection point** in a melatonin profile by calling
#' `seek_inflection()`, which performs a grid search over a **region of interest (ROI)**.
#'
#' @param profile_data A tibble containing melatonin profile data with columns:
#'   - `datetime` (POSIXct): The timestamp of each measurement.
#'   - `melatonin` (numeric): The melatonin concentration.
#'   - `base` (binary): Indicator for the baseline segment.
#'   - `ascending` (binary): Indicator for the ascending segment.
#'   - `intermediate` (optional, binary): Indicator for an intermediate segment.
#' @param threshold Numeric. The threshold value for constraining the fit (default = 2.3).
#' @param posix_roi A list defining the **region of interest** (ROI) using POSIXct timestamps:
#'   - `x_start`: Start time of the ROI.
#'   - `x_end`: End time of the ROI.
#'   - `y_min`: Minimum melatonin concentration in the ROI.
#'   - `y_max`: Maximum melatonin concentration in the ROI.
#' @param fit_type Character. Specifies the type of fit used:
#'   - `"linear"` (default): Uses a linear model.
#'   - `"parabolic"`: Uses a parabolic model for the ascending segment.
#' @param fine_flag Logical. If `TRUE`, enables **fine-grid refinement** in `seek_inflection()` (default = `TRUE`).
#'
#' @return A list containing:
#'   - `inflection_point`: The estimated inflection point (`x`, `y`).
#'   - `base_params`: Parameters of the **base segment** fit.
#'   - `ascending_params`: Parameters of the **ascending segment** fit.
#'   - `grid_big`: The coarse grid used in the first search.
#'   - `res_big`: The residuals corresponding to `grid_big`.
#'   - `grid_small`: The fine grid used in the second search (if `fine_flag = TRUE`).
#'   - `res_small`: The residuals corresponding to `grid_small` (if `fine_flag = TRUE`).
#'
#' @details
#' - The **ROI is converted** from POSIXct format to decimal time using `posixct_to_decimal()`.
#' - If the `profile_data` contains an `"intermediate"` segment, it is **included** in the fit.
#' - Calls `seek_inflection()` to identify the best-fit **inflection point**.
#'
#' @seealso
#' \code{\link{seek_inflection}}, \code{\link{posixct_to_decimal}}, \code{\link{fit}}
#'
#' @export
get_inflection <- function(profile_data, threshold = 2.3, posix_roi, fit_type = "linear", fine_flag = TRUE) {

  # Select relevant segments: include "intermediate" if it exists
  if ("intermediate" %in% colnames(profile_data)) {
    filtered_data <- dplyr::filter(profile_data, base == 1 | ascending == 1 | intermediate == 1)
  } else {
    filtered_data <- dplyr::filter(profile_data, base == 1 | ascending == 1)
  }

  if (nrow(filtered_data) < 3) {
    stop(
      "Cannot estimate DLMO because fewer than three points remain in the fitting region ",
      "after segment selection. This usually means the selected threshold/interval produced ",
      "too few base, intermediate, or ascending points. Inspect the segmented profile and consider ",
      "increasing `interval_limit` to include a later sustained rise, or adjusting `threshold` so ",
      "the selected rise better matches the intended DLMO event.",
      call. = FALSE
    )
  }

  datetime_ref <- filtered_data$datetime[3]
  if (is.na(datetime_ref)) {
    stop(
      "Cannot estimate DLMO because the fitting-region datetime reference is missing. ",
      "Inspect the segmented profile and consider adjusting `threshold` and/or `interval_limit`.",
      call. = FALSE
    )
  }

  # Convert POSIXct ROI coordinates to decimal time for numerical fitting.
  # Use the same datetime reference as fit(), which operates on filtered_data.
  roi <- list(
    x = posixct_to_decimal(c(posix_roi$x_start, posix_roi$x_end), datetime_ref),
    y = c(posix_roi$y_min, posix_roi$y_max)
  )
  attr(roi, "datetime_ref") <- datetime_ref

  if (!all(is.finite(c(roi$x, roi$y)))) {
    stop(
      "Cannot estimate DLMO because the search region could not be converted to finite numeric bounds. ",
      "This usually means the selected threshold/interval produced invalid base, intermediate, or ",
      "ascending segments. Inspect the segmented profile and consider increasing `interval_limit` ",
      "to include a later sustained rise, or adjusting `threshold` so the selected rise better ",
      "matches the intended DLMO event.",
      call. = FALSE
    )
  }

  # Run the inflection search using the selected data and ROI
  poi <- seek_inflection(filtered_data, threshold = threshold, roi, fit_type = fit_type, fine_flag = fine_flag)

  return(poi)
}
