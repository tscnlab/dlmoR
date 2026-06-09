#' Define Intermediate Segment of Melatonin Profile
#'
#' This function identifies and labels the intermediate segment of a melatonin profile.
#' The intermediate segment lies between the base segment and the ascending segment,
#' where melatonin concentrations begin to rise but have not crossed the threshold value.
#'
#' @param profile_data A tibble containing the melatonin profile with the following columns:
#'   - `datetime`: A POSIXct column representing timestamps for each measurement.
#'   - `melatonin`: Numeric column representing melatonin concentrations.
#'   - `base`: Binary column (1 for base segment, 0 otherwise).
#'   - `ascending`: Binary column (1 for ascending segment, 0 otherwise).
#' @param threshold Numeric. The melatonin concentration threshold for defining profile segments (default = 2.3 pg/mL).
#' @return A tibble with an additional `intermediate` column:
#'   - `intermediate`: Binary column (1 for intermediate segment, 0 otherwise).
#'
#' If no intermediate segment is identified, the `intermediate` column will not be created.
#'
#' @examples
#' library(dplyr)
#' library(lubridate)
#'
#' # Example data
#' profile_data <- tibble(
#'   datetime = seq(ymd_hms("2023-01-01 20:00:00"), by = "15 min", length.out = 12),
#'   melatonin = c(1.2, 1.4, 1.5, 2.0, 2.5, 2.3, 2.8, 3.5, 4.2, 4.7, 5.0, 5.2),
#'   base = c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0),
#'   ascending = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1)
#' )
#'
#' # Define the intermediate segment
#' updated_profile <- define_intermediate_segment(profile_data, threshold = 2.3)
#' print(updated_profile)
#' @export
define_intermediate_segment <- function(profile_data, threshold = 2.3) {
  # Ensure the data is sorted by datetime for consistent processing
  profile_data <- profile_data %>% dplyr::arrange(.data$datetime)

  # Find the last row of the base segment
  # last_base_row <- max(which(profile_data$base == 1), na.rm = TRUE)
  base_indices <- which(profile_data$base == 1)
  last_base_row <- if (length(base_indices) > 0) max(base_indices) else NA

  # Find the first row of the ascending segment
  # first_ascending_row <- min(which(profile_data$ascending == 1), na.rm = TRUE)
  ascending_indices <- which(profile_data$ascending == 1)
  first_ascending_row <- if (length(ascending_indices) > 0) min(ascending_indices) else NA

  # If no base is found and the first row has NA slope, mark it as base
  if (is.na(last_base_row) & is.na(profile_data$slope[1])) {
    profile_data$base[1] <- 1
    last_base_row <- 1
  }

  # Check if there are rows between the base and ascending segments
  if (last_base_row < first_ascending_row - 1) {
    # Identify rows that belong to the intermediate segment
    intermediate_rows <- (last_base_row + 1):(first_ascending_row - 1)

    # Add a new column 'intermediate', defaulting to 0
    profile_data <- profile_data %>%
      dplyr::mutate(intermediate = 0)

    # Mark the intermediate rows with 1
    profile_data$intermediate[intermediate_rows] <- 1
  } else {
    # If no intermediate segment exists, add the column with all 0s
    # profile_data <- profile_data %>%
    #   dplyr::mutate(intermediate = 0)
  }

  # Return the updated profile data with the new intermediate column
  return(profile_data)
}
