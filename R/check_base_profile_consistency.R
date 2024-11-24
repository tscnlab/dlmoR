#' Check Consistency of Base Segment in Melatonin Profile
#'
#' This function checks for consistency issues in the base segment of a melatonin profile.
#' It evaluates the following conditions:
#' - Presence of NA values in the base segment.
#' - Descents in melatonin concentration across the threshold of 2.3 ng/mL within the base segment.
#' - Large slope differences (greater than half the threshold value) in the base segment.
#'
#' The function prints warnings for any detected inconsistencies and returns a tibble
#' containing the full profile with additional columns to indicate where inconsistencies occurred.
#'
#' @param profile_data A tibble containing melatonin profile data with the following columns:
#'   - `datetime` (POSIXct): Timestamp for each data point.
#'   - `melatonin` (numeric): Melatonin concentrations.
#'   - `time` (numeric): Time in hours or another unit (optional, not used in this function).
#'   - `slope` (numeric): Slope of melatonin concentration between consecutive time points.
#'   - `base` (integer): Indicator for base segment (1 for base, 0 otherwise).
#'
#' @return A tibble containing the original profile data with three additional logical columns:
#'   - `warning_na`: TRUE if an NA value is present in the base segment.
#'   - `warning_threshold_descend`: TRUE if there is a descent across the threshold within the base segment.
#'   - `warning_large_diff`: TRUE if a slope exceeds half the threshold value within the base segment.
#'
#' @examples
#' # Example profile data
#' profile_data <- tibble::tibble(
#'   datetime = seq.POSIXt(from = as.POSIXct("2024-11-23 00:00:00"),
#'                         by = "hour", length.out = 10),
#'   melatonin = c(2.5, 2.2, 2.1, 2.3, 2.8, 2.9, 3.0, 2.7, 2.4, NA),
#'   time = seq(0, 9, 1),
#'   slope = c(-0.3, -0.1, 0.2, 0.5, 0.1, 0.1, -0.3, -0.3, NA, NA),
#'   base = c(1, 1, 1, 1, 1, 0, 0, 0, 0, 0)
#' )
#'
#' # Check for inconsistencies in the base segment
#' .check_base_profile_consistency(profile_data)
#'
#' @keywords internal
.check_base_profile_consistency <- function(profile_data) {
  # Check for NA values in the base segment
  warning_na <- dplyr::filter(profile_data, .data$base == 1) %>%
    dplyr::mutate(warning_na = is.na(.data$melatonin)) %>%
    dplyr::filter(warning_na)

  if (nrow(warning_na) > 0) {
    warning("Warning: NA values found in base segment at the following timestamps:\n",
            paste(warning_na$datetime, collapse = ", "))
  }

  # Check for descents across threshold in the base segment
  warning_threshold_descend <- dplyr::filter(profile_data, .data$base == 1) %>%
    dplyr::mutate(warning_threshold_descend = .data$melatonin > 2.3 & dplyr::lead(.data$melatonin, default = 0) <= 2.3) %>%
    dplyr::filter(warning_threshold_descend)

  if (nrow(warning_threshold_descend) > 0) {
    warning("Warning: Descents across threshold found in base segment at the following timestamps:\n",
            paste(warning_threshold_descend$datetime, collapse = ", "))
  }

  # Check for large slope differences in the base segment
  warning_large_diff <- dplyr::filter(profile_data, .data$base == 1) %>%
    dplyr::mutate(warning_large_diff = abs(.data$slope) > 0.5 * 2.3) %>%
    dplyr::filter(warning_large_diff)

  if (nrow(warning_large_diff) > 0) {
    warning("Warning: Large slope differences found in base segment at the following timestamps:\n",
            paste(warning_large_diff$datetime, collapse = ", "))
  }

  # Combine all warnings into a single tibble
  warnings_combined <- dplyr::left_join(
    profile_data,
    dplyr::bind_rows(
      warning_na %>% dplyr::select(.data$datetime, warning_na),
      warning_threshold_descend %>% dplyr::select(.data$datetime, warning_threshold_descend),
      warning_large_diff %>% dplyr::select(.data$datetime, warning_large_diff)
    ),
    by = "datetime"
  ) %>%
    dplyr::mutate(
      warning_na = dplyr::coalesce(warning_na, FALSE),
      warning_threshold_descend = dplyr::coalesce(warning_threshold_descend, FALSE),
      warning_large_diff = dplyr::coalesce(warning_large_diff, FALSE)
    )

  return(warnings_combined)
}

