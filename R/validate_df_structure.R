#' Function to Validate Data
#'
#' This function validates the structure of an input dataframe. It verifies that the dataframe
#'     contains `melatonin` and `datetime` (or `time`) columns and checks that they contain numeric and POSIXct (or hms) data, respectively.
#'
#' @param data A data frame, ideally with columns `datetime` (or `time`) and `melatonin`.
#' @return Returns a dataframe with columns `melatonin`, `datetime` and `time` if the input df contains
#'    `melatonin` and `datetime` columns. Returns `melatonin` and `time` columns if the input df only contains these.
#'    If criteria are not met, a corresponding error message is returned.
#'
validate_df_structure <- function(data) {
  # Ensure the data is in tibble format
  data <- tidyr::as_tibble(data)
  # Required columns
  required_columns <- c("melatonin")

  # Check for required columns
  if (!all(required_columns %in% colnames(data))) {
    stop("The tibble must contain a 'melatonin' column.")
  }

  # Handle 'time' column with POSIXct data
  if ("time" %in% colnames(data) && inherits(data$time, "POSIXct")) {
    message("'time' column contains POSIXct data; renaming to 'datetime' and extracting time component into a new 'time' column.")
    data <- data %>%
      dplyr::rename(datetime = .data$time) %>%
      dplyr::mutate(time = hms::as_hms(.data$datetime))
  }

  # Handle 'datetime' column with hms data
  if ("datetime" %in% colnames(data) && inherits(data$datetime, "hms")) {
    colnames(data)[colnames(data) == "datetime"] <- "time"
    stop("The 'datetime' column contains hms data. It has been renamed to 'time'.")
  }

  # Check for 'datetime' column and validate POSIXct
  if ("datetime" %in% colnames(data)) {
    if (!inherits(data$datetime, "POSIXct")) {
      stop("The 'datetime' column must be of class 'POSIXct'.")
    }
    # Create 'time' column if not present
    if (!"time" %in% colnames(data)) {
      message("Creating 'time' column from 'datetime'.")
      data <- data %>% dplyr::mutate(time = hms::as_hms(.data$datetime))
      #print(data)
    }
  }

  # Check for 'time' column and validate hms format
  if ("time" %in% colnames(data)) {
    if (!inherits(data$time, "hms")) {
      stop("The 'time' column must be of class 'hms'. If using a character column, convert it to 'hms' or 'POSIXct' format first.")
    }
  } else if (!"datetime" %in% colnames(data)) {
    stop("The tibble must contain either a 'time' column of class 'hms' or a 'datetime' column of class 'POSIXct'.")
  }

  # Ensure 'melatonin' column is numeric
  if (!is.numeric(data$melatonin)) {
    stop("The 'melatonin' column must contain numeric values.")
  }

  return(data)
}
