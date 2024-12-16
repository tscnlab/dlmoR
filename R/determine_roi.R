# define region of interest
define_roi <- function(profile_data, threshold = 2.3) {

  # Ensure the profile_data is sorted by datetime
  profile_data <- profile_data %>% dplyr::arrange(datetime)

  # Extract relevant points based on base and ascending flags
  base_points <- profile_data %>% dplyr::filter(base == 1)
  ascending_points <- profile_data %>% dplyr::filter(ascending == 1)
  intermediate_points <- profile_data %>% dplyr::filter(intermediate == 1)
  print("intermediate"%in%colnames(profile_data))
  if ("intermediate"%in%colnames(profile_data)){
    x_start <- intermediate_points$datetime[1]
    print("x-start")
    print(x_start)
  } else if (nrow(base_points) < 2) {
    # If the base segment is a single node
    x_start <- base_points$datetime[1]
  } else {
    # Otherwise, take the midpoint of the last two base points
    x_start<- decimal_to_posixct(mean(posixct_to_decimal(base_points$datetime[(nrow(base_points) - 1):nrow(base_points)], profile_data$datetime)), profile_data$datetime)
    #x_start <- as.POSIXct(mean(as.numeric(base_points$datetime[(nrow(base_points) - 1):nrow(base_points)])), origin = "1970-01-01")
  }

  if (nrow(profile_data %>% dplyr::filter(base == 0 & ascending == 0)) > 0) { #TODO FIX THIS TO INCORPORATE INTERMEDIATE
    # Include intermediate segments and up to 95% of the first ascending segment
    # x_end <- ascending_points$datetime[1] + 0.95 * (ascending_points$datetime[2] - ascending_points$datetime[1])
    x_end <- tail(base_points$datetime,n=1) + 0.95 * (ascending_points$datetime[1] - tail(base_points$datetime,n=1))
  } else {
    # If no intermediate segments
    # x_end <- ascending_points$datetime[1] + 0.95 * (ascending_points$datetime[2] - ascending_points$datetime[1])
    x_end <- tail(base_points$datetime,n=1) + 0.95 * (ascending_points$datetime[1] - tail(base_points$datetime,n=1))

  }

  # Define the vertical bounds
  y_min <- min(profile_data$melatonin[profile_data$datetime >= x_start & profile_data$datetime <= x_end])
  y_max <- threshold

  list(
    x_start = x_start,
    x_end = x_end,
    y_min = y_min,
    y_max = y_max
  )
}
