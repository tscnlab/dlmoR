#' #' Title
#' #' @param profile a vector of melatonin concentrations
#' #' @return the segment of the melatonin concentrations defined as belonging to the baseline
#' #' @export
#' #'
#' #' @examples
#' # Function to define the base segment
#' # define_base_segment <- function(profile, threshold=2.3) {
#' #   # Calculate slopes
#' #   slopes <- calculate_slopes(profile)
#' #
#' #   # Create a tibble for the profile with slopes and other computations
#' #   profile_data <- tibble(
#' #     index = seq_along(profile),                # Index for reference
#' #     melatonin = profile,                       # Melatonin concentration values
#' #     slope = c(slopes, NA),                     # Append NA to align slope with profile length
#' #     base = NA                                  # Initialize base column
#' #   )
#' #
#' #   # Identify base segment based on slope and threshold
#' #   profile_data <- profile_data %>%
#' #     mutate(
#' #       base = if_else(
#' #         slope <= 0 & (melatonin <= threshold | lead(melatonin, default = NA) <= threshold),
#' #         1,
#' #         0,
#' #         missing = 0
#' #       )
#' #     ) %>%
#' #     # Propagate the base segment until the rightmost identified base
#' #     mutate(
#' #       base = if_else(row_number() <= max(which(base == 1), na.rm = TRUE) + 1, 1, base, missing = base)
#' #     )
#' #
#' #   # Filter the base segment and return it
#' #   base_segment <- profile_data %>% filter(base == 1)
#' #
#' #   # Return the tibble containing the identified base segment
#' #   return(base_segment)
#' # }
#' define_base_segment <- function(profile, threshold = 2.3) {
#'   # Calculate slopes
#'   slopes <- calculate_slopes(profile)
#'
#'   # Create a tibble for the full profile with slopes and base labels
#'   profile_data <- tibble(
#'     index = seq_along(profile),                # Index for x-axis
#'     melatonin = profile,                       # Melatonin concentration values
#'     slope = c(slopes, NA),                     # Append NA to align slope with profile length
#'     base = rep(0, length(profile))             # Initialize base column to 0 for all points
#'   )
#'
#'   # Identify base segment based on slope and threshold
#'   profile_data <- profile_data %>%
#'     mutate(
#'       base = if_else(
#'         slope <= 0 & (melatonin <= threshold | lead(melatonin, default = NA) <= threshold),
#'         1,
#'         base,
#'         missing = base
#'       )
#'     ) %>%
#'     # Propagate the base segment until the rightmost identified base
#'     mutate(
#'       base = if_else(row_number() <= max(which(base == 1), na.rm = TRUE) + 1, 1, base, missing = base)
#'     )
#'
#'   # Return the full tibble with base labels
#'   return(profile_data)
#' }
#'
#' #' Define the base segment of melatonin profile
#' @param profile_tibble A tibble containing datetime, melatonin, and time columns
#' @param threshold The melatonin concentration threshold to define the base segment (default = 2.3)
#' @return A tibble with the base segment labeled (1 for base, 0 for non-base), including datetime, time, and melatonin columns
#' @export
#'
#' @examples
define_base_segment <- function(profile_tibble, threshold = 2.3) {
  # Extract melatonin values from the input tibble (profile_tibble)
  profile <- profile_tibble$melatonin

  # Calculate slopes based on the melatonin profile
  slopes <- calculate_slopes(profile)

  # Create a tibble for the full profile with slopes and base labels
  profile_data <- profile_tibble %>%
    dplyr::mutate(
      slope = c(slopes, NA),                     # Append NA to align slope with profile length
      base = rep(0, length(profile))              # Initialize base column to 0 for all points
    )

  # Identify base segment based on slope and threshold
  profile_data <- profile_data %>%
    dplyr::mutate(
      base = dplyr::if_else(
        .data$slope <= 0 & (.data$melatonin <= threshold | dplyr::lead(.data$melatonin, default = NA) <= threshold),
        1,  # Mark base as 1
        .data$base,  # Keep current base value
        missing = .data$base  # Keep NA for missing values
      )
    ) %>%
    # Propagate the base segment until the rightmost identified base
    dplyr::mutate(
      base = dplyr::if_else(dplyr::row_number() <= max(which(.data$base == 1), na.rm = TRUE) + 1, 1, .data$base, missing = .data$base)
    )

  # Return the full tibble with datetime, melatonin, time, slope, and base columns
  return(profile_data)
}

