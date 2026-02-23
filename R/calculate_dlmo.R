#' Calculate Dim-Light Melatonin Onset (DLMO)
#'
#' This function calculates the DLMO time point based on input melatonin concentration data
#' with an associated time series. The user can provide data directly as a data frame
#' or specify a .csv file to load the data.
#'
#' @param data A data frame with columns `datetime` (POSIXct) and `melatonin` (numeric). If NULL, use
#' the `file_path` parameter to load data.
#' @param file_path A string specifying the path to a CSV file containing the data.
#' The file must have two columns: `datetime` (POSIXct) and `melatonin` (numeric).
#' @param threshold The numeric melatonin threshold for defining DLMO (default: 2.3 pg/mL).
#' @param interval_limit Numeric or lubridate duration. Indicates the time window within which,
#' if two melatonin threshold crossings occur, the second crossing is taken to represent
#' the melatonin rise. If provided as a numeric value, it represents the window in hours
#' (e.g., `interval_limit = 2` for 2 hours, or `interval_limit = 0.5` for 30 minutes).
#' Alternatively, users can specify a `lubridate` duration object
#' (e.g., `lubridate::hours(2)` for 2 hours, `lubridate::minutes(30)` for 30 minutes).
#' Default is 2 hours.
#' @param fine_flag Logical. If `TRUE`, performs an additional fine-grid search to refine the DLMO point after the initial coarse search (default: `TRUE`).
#' @return A list containing the following elements:
#'
#' - **`prof`**: A tibble containing the segmented melatonin profile.
#'   - `datetime` (POSIXct): Timestamps of melatonin measurements.
#'   - `melatonin` (numeric): Melatonin concentration values.
#'   - `time` (hms): Time of day.
#'   - `slope` (numeric): Rate of change in melatonin concentration.
#'   - `base` (binary, 0/1): Indicates baseline segment.
#'   - `ascending` (binary, 0/1): Indicates ascending segment.
#'   - `intermediate` (binary, 0/1): if present, indicates intermediate segment.
#'
#' - **`prl`**: List containing parallelogram rule fit paramters.
#'   - `pll_datetime_0` (POSIXct): Estimated lower bound of melatonin rise.
#'   - `pll_datetime_1` (POSIXct): Estimated upper bound of melatonin rise.
#'   - `pll_slope` (numeric): Slope of melatonin rise.
#'   - `corners` (list): Coordinates of parallelogram corners:
#'     - `ll`, `lr`, `ur`, `ul` (numeric vectors): (x, y) for each corner.
#'   - `flag` (logical): Whether parallelogram rule was violated or not.
#'
#' - **`roi`**: List defining the Region of Interest (ROI).
#'   - `x_start` (POSIXct): Start time of the ROI window.
#'   - `x_end` (POSIXct): End time of the ROI window.
#'   - `y_min` (numeric): Minimum melatonin value in the ROI.
#'   - `y_max` (numeric): Maximum melatonin value in the ROI.
#'
#' - **`ip`**: List containing estimated DLMO inflection point.
#'   - `inflection_point` (tibble):
#'     - `x` (numeric): Estimated DLMO time index in units of decimal-hours.
#'     - `y` (numeric): Melatonin concentration at DLMO.
#'   - `base_params` (numeric): Parameter for linear-fit of base segment of profile.
#'   - `ascending_params` (list): Parameters for the linear or parabolic fit of melatonin rise:
#'     - `a`, `b`, `c` (numeric): Fitted parameters.
#' - **`dlmo_time`**: Estimated DLMO time-stamp in units of hh:mm:ss
#' - **`dlmoplotcoarse`**: `ggplot` object visualizing melatonin profile (coarse view).
#' - **`dlmoplotfine`**: `ggplot` object visualizing melatonin profile (fine view).
#'
#' @details
#' **Understanding the Outputs:**
#'
#' - The `prof` tibble contains labeled data for different melatonin profile phases.
#' - The `prl` list provides the parallelogram rule fit used to trim the melatonin rise segment of the profile to ensure only strong rises are fit when determining DLMO.
#' - The `roi` specifies the bounds of the search window used for DLMO determination.
#' - The `ip` list contains the  DLMO inflection point, the fit parameters for the base and ascending regions, as well as the grid of residuals from fitting the melatonin profile at each point of the ROI.
#' - The `dlmo_time` variable contains the DLMO timestamp in units of hh:mm:ss (equivalent to `ip$inflection_point$x` in decimal-hours)
#' - The plots (`dlmoplotcoarse` and `dlmoplotfine`) visualize the results of the coarse and fine grid DLMO search, respectively.


#' @examples
#' # Load the sample melatonin profile data included in the package
#' filename <- system.file("extdata/sample_melatonin_profile.csv", package = "dlmoR")
#'
#' # Calculate the DLMO using the sample data and a threshold of 5
#' sample_dlmo <- calculate_dlmo(file_path = filename, threshold = 5)
#'
#' @export
#'
calculate_dlmo <- function(data = NULL, file_path = NULL, threshold = 2.3, interval_limit = lubridate::hours(2), fine_flag = TRUE) {

  # Check if input is provided either directly or via file
  if (is.null(data) && is.null(file_path)) {
    stop("You must provide either `data` or `file_path`.")
  }

  # Convert numeric interval_limit to a lubridate duration (interpreted as hours)
  if (is.numeric(interval_limit)) {
    interval_limit <- lubridate::hours(interval_limit)
  } else if (!lubridate::is.duration(interval_limit)) {
    stop("`interval_limit` must be a numeric value (interpreted as hours) or a lubridate duration.")
  }
  # Load data if file_path is provided
  if (!is.null(file_path)) {
    message("Loading data from file: ", file_path)
    data <- .read_melatonin_data(file_path)
  }
  # Validate df structure
  data<-validate_df_structure(data)


  # Validate melatonin profile
  prf<-preprocess_profile(data, threshold = threshold)

  # define & truncate profile segments
  prf<-define_base_segment(prf, threshold = threshold)
  .check_base_profile_consistency(prf, threshold = threshold)
  prf<-define_ascending_segment(prf, threshold = threshold, interval_limit = interval_limit)
  prf<-truncate_ascending_segment(prf)
  prf$profile<-truncate_base_segment(prf$profile, threshold = threshold)
  prf$profile<-define_intermediate_segment(prf$profile, threshold = threshold)

  # define roi & search for dlmo inflection point
  roix<-define_roi(profile_data = prf$profile, threshold = threshold)
  ipx <- get_inflection(prf$profile, threshold = threshold, roix, fine_flag = fine_flag)

  # Initialize as NULL
  ipx_coarse <- NULL
  ipx_fine <- NULL

  # Assign only the relevant one based on fine_flag
  if (fine_flag) {
    ipx_fine <- list(
      inflection_point = ipx$inflection_point_fine,
      base_params = ipx$base_params_fine,
      ascending_params = ipx$ascending_params_fine,
      grid = ipx$grid_small,
      res = ipx$res_small
    )

    ipx_coarse <- list(
      inflection_point = ipx$inflection_point_coarse,
      base_params = ipx$base_params_coarse,
      ascending_params = ipx$ascending_params_coarse,
      grid = ipx$grid_big,
      res = ipx$res_big
    )
  } else {
    ipx_coarse <- list(
      inflection_point = ipx$inflection_point_coarse,
      base_params = ipx$base_params_coarse,
      ascending_params = ipx$ascending_params_coarse,
      grid = ipx$grid_big,
      res = ipx$res_big
    )
  }



  # Initialize dlmo list
  dlmo <- list(coarse = NULL, fine = NULL)

  # Process coarse if available
  if (!is.null(ipx_coarse)) {
    dlmo_coarse_time <- hms::as_hms(decimal_to_posixct(ipx_coarse$inflection_point$x, prf$profile$datetime))
    dlmo$coarse <- list(
      time = dlmo_coarse_time,
      fit_melatonin = ipx_coarse$inflection_point$y
    )

    # cat("\nCOARSE ascending_params:\n")
    # str(ipx_coarse$ascending_params)
    #
    # cat("\nFINE ascending_params:\n")
    # if (!is.null(ipx_fine)) str(ipx_fine$ascending_params)

    dlmo$coarse <- process_fits(ipx_coarse, dlmo$coarse)
  }

  # Process fine if available
  if (!is.null(ipx_fine)) {
    dlmo_fine_time <- hms::as_hms(decimal_to_posixct(ipx_fine$inflection_point$x, prf$profile$datetime))
    dlmo$fine <- list(
      time = dlmo_fine_time,
      fit_melatonin = ipx_fine$inflection_point$y
    )
    dlmo$fine <- process_fits(ipx_fine, dlmo$fine)
  }


  # create and save visualizations
  vis_coarse<-plot_profile(prf$profile, show_threshold = TRUE, threshold = threshold, show_segments = TRUE, show_parallelogram = TRUE, pll_result = prf$plll, show_roi = TRUE, roi = roix, show_dlmoIP = TRUE, dlmoFit = ipx_coarse, dlmo = dlmo, show_fit = TRUE, show_roi_heatmap = TRUE, plot_coarse = TRUE)
  if(fine_flag){
  vis_fine<-plot_profile(prf$profile, show_threshold = TRUE, threshold = threshold, show_segments = TRUE, show_parallelogram = TRUE, pll_result = prf$plll, show_roi = TRUE, roi = roix, show_dlmoIP = TRUE, dlmoFit = ipx_fine, dlmo = dlmo, show_fit = TRUE, show_roi_heatmap = TRUE, plot_coarse = FALSE)
  }
  else{
    vis_fine <- NULL
  }
  return(list(prof = prf$profile, prl = prf$plll, roi = roix, ip = ipx, dlmo = dlmo, dlmoplotcoarse = vis_coarse, dlmoplotfine = vis_fine))
}

#' Helper Function to Read-in Melatonin Data from a CSV-File
#'
#' This function reads in a CSV file which ideally contains melatonin concentration data and associated datetime (or time) stamps and returns
#' a tibble.
#'
#' @param file_path A string specifying the path to the CSV file.
#' @return A data frame with columns `datetime` (or `time` and `melatonin`).
#' @keywords internal
.read_melatonin_data <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("The file does not exist: ", file_path)
  }
  # data <- read.csv(file_path, sep = ";", header = TRUE)
  data <- readr::read_delim(file_path, delim = NULL, show_col_types = FALSE)
  # # check if file is a dataframe
  # if (!is.data.frame(data)) {
  #   stop("`data` must be a data frame containing `datetime` (or `time`) and `melatonin` columns.")
  # }

  return(data)
}

# Internal function to calculate DLMO
.find_dlmo <- function(time, melatonin, threshold) {
  # Find the first time point where melatonin exceeds the threshold
  index <- which(melatonin >= threshold)[1]
  if (is.na(index)) {
    stop("No DLMO detected: melatonin does not exceed the threshold.")
  }
  return(time[index])
}

process_fits <- function(ip, dlmo) {
  # Initialize fit_lines list
  dlmo$fit_lines <- list()

  # Ensure inflection point exists and has x & y
  if (!is.null(ip$inflection_point) &&
      !is.null(ip$inflection_point$x) &&
      !is.null(ip$inflection_point$y)) {

    # Extract x and y from inflection point
    inf_x <- ip$inflection_point$x
    inf_y <- ip$inflection_point$y

    # Compute baseline fit intercept
    base_intercept <- inf_y - ip$base_params * inf_x

    # Store base fit with computed intercept
    dlmo$fit_lines$base <- list(
      type = "linear",
      m = ip$base_params,   # Slope
      b = base_intercept    # Computed intercept
    )

    # Check if ascending fit is linear or parabolic
    if (length(ip$ascending_params) == 1) {
      # Compute ascending intercept
      ascending_intercept <- inf_y - ip$ascending_params * inf_x

      dlmo$fit_lines$ascending <- list(
        type = "linear",
        m = ip$ascending_params,  # Slope
        b = ascending_intercept   # Computed intercept
      )
    } else if (length(ip$ascending_params) == 3) {
      # Store parabolic fit (intercept already given as 'c')
      dlmo$fit_lines$ascending <- list(
        type = "parabolic",
        a = ip$ascending_params$a,
        b = ip$ascending_params$b,
        c = ip$ascending_params$c
      )
    }

  } else {
    # Handle missing inflection point
    warning("Missing or incomplete ip$inflection_point! Storing base fit without intercept.")

    # Store base without computing intercept
    dlmo$fit_lines$base <- list(
      type = "linear",
      m = ip$base_params,   # Store only slope if intercept can't be calculated
      b = NA                # Set b as NA since we couldn't compute it
    )

    # Check if ascending fit is linear or parabolic
    if (length(ip$ascending_params) == 1) {
      dlmo$fit_lines$ascending <- list(
        type = "linear",
        m = ip$ascending_params,
        b = NA  # Can't compute without inflection point
      )
    } else if (length(ip$ascending_params) == 3) {
      dlmo$fit_lines$ascending <- list(
        type = "parabolic",
        a = ip$ascending_params$a,
        b = ip$ascending_params$b,
        c = ip$ascending_params$c
      )
    }
  }

  # Return the updated dlmo object
  return(dlmo)
}

