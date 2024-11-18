read_melatonin_data <- function(file_path) {
  data <- read.csv(file_path)
  # Validate required columns
  if (!all(c("time", "melatonin") %in% colnames(data))) {
    stop("CSV file must contain 'time' and 'melatonin' columns.")
  }
  return(data)
}
