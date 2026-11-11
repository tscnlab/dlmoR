#' Truncate the Base Segment of a Melatonin Profile
#'
#' This function identifies and removes points in the base segment where melatonin
#' concentrations exceed a specified threshold. It ensures that only the earliest
#' base points remain while maintaining a valid base segment structure.
#'
#' @param profile_data A tibble or dataframe containing the melatonin profile.
#'   It must include the columns:
#'   \itemize{
#'     \item `datetime` (POSIXct) – Timestamps of data points.
#'     \item `melatonin` (numeric) – Melatonin concentration levels.
#'     \item `base` (integer, 0 or 1) – Indicator of base segment points.
#'   }
#' @param threshold Numeric. The melatonin concentration threshold. Any base segment
#'   points above this threshold are removed. Default is `2.3`.
#'
#' @return A modified tibble where invalid base points have been removed.
#'   If no valid base points remain, a warning is issued.
#' @export
truncate_base_segment <- function(profile_data, threshold = 2.3) {

  # Ensure the data is sorted by datetime before processing
  profile_data <- profile_data %>% dplyr::arrange(.data$datetime)

  # Identify and remove base segment points above the threshold at the start of the base segment
  profile_data <- profile_data %>%
    dplyr::mutate(
      base = dplyr::if_else(
        .data$base == 1 &
          # Ensure `NA` values are handled properly before `cumsum()`
          cumsum(dplyr::coalesce(.data$base == 1 & .data$melatonin <= threshold, FALSE)) == 0,
        0,  # Set base to 0 for these invalid points
        .data$base
      )
    )

  # Ensure valid base points exist before checking `all(profile_data$base == 0)`
  if (nrow(profile_data) == 0 || all(is.na(profile_data$base))) {
    warning("No valid data found in profile_data.")
    return(profile_data)
  }

  # Check if all base points have been removed
  if (all(profile_data$base == 0)) {
    warning("Hockey-stick time: no base part")  # Notify that no base segment remains
    return(profile_data)  # Return profile data unchanged
  }

  # Return the updated profile data with truncated base points
  return(profile_data)
}
