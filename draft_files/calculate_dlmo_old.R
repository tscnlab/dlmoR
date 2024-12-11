#' Calculate Dim-Light Melatonin Onset (DLMO)
#'
#' This function calculates the DLMO based on input melatonin concentration data
#' with an associated time series. The user can provide data directly as a data frame
#' or specify a file to load the data.
#'
#' @param data A data frame with columns `time` and `melatonin`. If NULL, use
#' the `file_path` parameter to load data.
#' @param file_path A string specifying the path to a CSV file containing the data.
#' The file must have two columns: `time` (numeric or POSIXct) and `melatonin` (numeric).
#' @param threshold Numeric. The melatonin threshold for defining DLMO (default: 10).
#' @return A numeric or time value indicating the calculated DLMO time.
#' @export
#'
#' @examples
#' # Using a data frame directly:
#' df <- data.frame(time = seq(1, 6), melatonin = c(5, 7, 9, 11, 13, 15))
#' calculate_dlmo(df)
#'
#' # Loading data from a file:
#' filename<-system.file("extdata/CiViBe_204_FD.csv", package = "dlmoR")
#' # calculate_dlmo(file_path = "melatonin_data.csv")
#' @export
calculate_dlmo <- function(data = NULL, file_path = NULL, threshold = 2.3) {
  # Check if input is provided either directly or via file
  if (is.null(data) && is.null(file_path)) {
    stop("You must provide either `data` or `file_path`.")
  }

  # Load data if file_path is provided
  if (!is.null(file_path)) {
    message("Loading data from file: ", file_path)
    data <- read_melatonin_data(file_path)
  }

  # Validate data structure
  data<-validate_df_structure(data)
  # print(data1)
  return(data)

  # # Validate and process the data
  # if (!is.data.frame(data)) {
  #   stop("`data` must be a data frame containing `datetime` (or `time`) and `melatonin` columns.")
  # }

  # if (!all(c("time", "melatonin") %in% colnames(data))) {
  #   stop("The data frame must contain columns `time` and `melatonin`.")
  # }
  #
  # if (!is.numeric(data$melatonin)) {
  #   stop("The `melatonin` column must contain numeric values.")
  # }
  #
  # # Check for valid time series
  # if (!is.numeric(data$time) && !inherits(data$time, "POSIXct")) {
  #   stop("The `time` column must be numeric or a POSIXct time series.")
  # }

  # # Perform the DLMO calculation
  # dlmo_time <- .find_dlmo(data$time, data$melatonin, threshold)
  # return(dlmo_time)
}

#' Helper Function to Read and Validate Melatonin Data from a File
#'
#' This function reads a CSV file containing melatonin concentration data and validates
#' its structure. The file must contain `time` and `melatonin` columns.
#'
#' @param file_path A string specifying the path to the CSV file.
#' @return A validated data frame with columns `time` and `melatonin`.
#' @keywords internal
read_melatonin_data <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("The file does not exist: ", file_path)
  }
  # data <- read.csv(file_path, sep = ";", header = TRUE)
  data <- read_delim(file_path, delim = NULL)
  # print(data)
  # if (!all(c("time", "melatonin") %in% colnames(data))) {
  #   stop("The file must contain columns 'time' and 'melatonin'.")
  # }
  # if (!is.numeric(data$melatonin)) {
  #   stop("The `melatonin` column in the file must contain numeric values.")
  # }
  # if (!is.numeric(data$time) && !inherits(data$time, "POSIXct")) {
  #   stop("The `time` column in the file must be numeric or a POSIXct time series.")
  # }

  # check if file is a dataframe
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame containing `datetime` (or `time`) and `melatonin` columns.")
  }

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

#' Function to Validate Data
#'
#' This function validates the structure of an input dataframe. It verifies that the dataframe
#'    it contains `melatonin` and `time` columns and checks that they contain numeric and/or POSIXct data.
#'
#' @param data A data frame with columns `time` and `melatonin`.
#' @return Returns a Boolean. If criteria are not met, an error message is returned.
validate_df_structure <- function(data) {
  # if (!all(c("time", "melatonin") %in% colnames(data))) {
  #   stop("The file must contain columns 'time' and 'melatonin'.")
  # }
  # if (!is.numeric(data$melatonin)) {
  #   stop("The `melatonin` column in the file must contain numeric values.")
  # }
  # if (!is.numeric(data$time) && !inherits(data$time, "POSIXct")) {
  #   stop("The `time` column in the file must be numeric or a POSIXct time series.")
  # }
  # Ensure the data is in tibble format
  data <- as_tibble(data)

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
      rename(datetime = time) %>%
      mutate(time = as_hms(datetime))
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
      data <- data %>% mutate(time = as_hms(datetime))
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

library(dplyr)
library(hms)

process_tibble <- function(data) {
  # Ensure the data is in tibble format
  data <- as_tibble(data)

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
      rename(datetime = time) %>%
      mutate(time = as_hms(datetime))
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
      data <- data %>% mutate(time = as_hms(datetime))
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

