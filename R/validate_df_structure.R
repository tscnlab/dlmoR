#' Validate and Standardize Dataframe Structure for Melatonin Analysis
#'
#' This function checks whether the input dataframe contains the required `melatonin`
#' column and at least one valid time reference (`datetime` or `time`). It ensures that:
#' - `melatonin` is numeric.
#' - `datetime` (if present) is in `POSIXct` format.
#' - `time` (if present) is in `hms` format.
#' - If only `time` is provided (as `POSIXct`), it is renamed to `datetime`, and a `time`
#'   column is created from it.
#' - If `datetime` exists but is of class `hms`, it is renamed to `time`, and an error is raised.
#'
#' @param data A dataframe (or tibble) expected to contain at least:
#'   \itemize{
#'     \item `melatonin` (numeric) – Melatonin concentration levels.
#'     \item `datetime` (POSIXct) – Timestamps, or
#'     \item `time` (hms or POSIXct) – Time of measurement.
#'   }
#'
#' @return A modified dataframe with:
#'   \itemize{
#'     \item `melatonin` (numeric)
#'     \item `datetime` (POSIXct) if available or converted
#'     \item `time` (hms) if `datetime` was present and converted
#'   }
#'
#' @details
#' If `time` is in `POSIXct`, it is renamed to `datetime`, and `time` is extracted from it.
#' If `datetime` is of class `hms`, it is renamed to `time`, and the function stops execution.
#' If neither `datetime` nor `time` exist, an error is raised.
#'
#' @examples
#' \dontrun{
#' df <- tibble::tibble(
#'   datetime = as.POSIXct(c("2024-04-16 12:00:00", "2024-04-16 12:30:00")),
#'   melatonin = c(1.2, 1.5)
#' )
#' validate_df_structure(df)
#' }
#' @export
validate_df_structure <- function(data) {
  # Convert input data to tibble format for consistency
  data <- tidyr::as_tibble(data)

  # Required column: 'melatonin'
  if (!"melatonin" %in% colnames(data)) {
    stop("The tibble must contain a 'melatonin' column.")
  }

  # Handle case where 'time' exists as POSIXct
  if ("time" %in% colnames(data) && inherits(data$time, "POSIXct")) {
    message("'time' column detected as POSIXct. Renaming to 'datetime' and extracting time.")
    data <- data %>%
      dplyr::rename(datetime = .data$time) %>%
      dplyr::mutate(time = hms::as_hms(.data$datetime))
  }

  # Handle case where 'datetime' exists as hms
  if ("datetime" %in% colnames(data) && inherits(data$datetime, "hms")) {
    colnames(data)[colnames(data) == "datetime"] <- "time"
    stop("The 'datetime' column is of class 'hms' and has been renamed to 'time'. Stopping execution.")
  }

  # Ensure 'datetime' column is in POSIXct format
  if ("datetime" %in% colnames(data)) {
    if (!inherits(data$datetime, "POSIXct")) {
      stop("The 'datetime' column must be of class 'POSIXct'.")
    }
    # If 'time' is missing, extract it from 'datetime'
    if (!"time" %in% colnames(data)) {
      message("Creating 'time' column from 'datetime'.")
      data <- data %>% dplyr::mutate(time = hms::as_hms(.data$datetime))
    }
  }

  # Ensure 'time' column (if present) is in hms format
  if ("time" %in% colnames(data)) {
    if (!inherits(data$time, "hms")) {
      stop("The 'time' column must be of class 'hms'. Convert it before running this function.")
    }
  } else if (!"datetime" %in% colnames(data)) {
    stop("The tibble must contain either a 'time' column (hms) or a 'datetime' column (POSIXct).")
  }

  # Ensure 'melatonin' column is numeric
  if (!is.numeric(data$melatonin)) {
    stop("The 'melatonin' column must contain numeric values.")
  }

  return(data)
}
