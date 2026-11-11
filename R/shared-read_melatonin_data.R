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
