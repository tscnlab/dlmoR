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
calculate_dlmo <- function(data = NULL, file_path = NULL, threshold = 10) {
  # Check if input is provided either directly or via file
  if (is.null(data) && is.null(file_path)) {
    stop("You must provide either `data` or `file_path`.")
  }

  # Load data if file_path is provided
  if (!is.null(file_path)) {
    message("Loading data from file: ", file_path)
    data <- .read_melatonin_data(file_path)
  }

  # Validate data structure
  data<-validate_df_structure(data)
  # print(data1)
  return(data)
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
  data <- read_delim(file_path, delim = NULL)

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



# library(dplyr)
# library(hms)


