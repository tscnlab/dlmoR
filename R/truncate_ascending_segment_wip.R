truncate_ascending_segment <- function(profile_data) {
  #print("allo")
  # Ensure there is an ascending segment to work with
  ascending_data <- profile_data %>% dplyr::filter(.data$ascending == 1)

  # If no ascending segment exists, issue a warning and return the data unchanged
  if (nrow(ascending_data) == 0) {
    warning("No ascending segments found in the profile data.")
    return(profile_data)  # Return unchanged if no ascending segment
  }


  # Find the steepest slope in the ascending segment (ignoring NA values)
  max_slope <- max(profile_data$slope[profile_data$ascending == 1], na.rm = TRUE)

  # Iteratively truncate until all rules are satisfied
  while (!check_rules(profile_data, max_slope)) {
    #print(profile_data, n = 22)
    # Identify the last ascending index
    last_ascending_index <- max(which(profile_data$ascending == 1))

    # Set the last point in the ascending segment to not ascending
    profile_data <- profile_data %>%
      dplyr::mutate(
        ascending = dplyr::if_else(
          dplyr::row_number() == last_ascending_index,
          0,
          .data$ascending
        )
      )
  }
  # print("before rule 3")
  # print(profile_data, n =28)
  # Find the steepest slope in the ascending segment (ignoring NA values)
  max_slope <- max(profile_data$slope[profile_data$ascending == 1], na.rm = TRUE)

  # Extract the two rightmost points in the ascending segment
  # Identify the last ascending index
  last_ascending_index <- max(which(profile_data$ascending == 1))
  #second_last_ascending_index<- last_ascending_index - 1

  # rightmost_points <- ascending_data %>%
  #   dplyr::arrange(desc(dplyr::row_number())) %>%
  #   head(2)
  #
  # # Calculate the slope of the rightmost segment
  # rightmost_slope <- (rightmost_points$melatonin[2] - rightmost_points$melatonin[1]) /
  #   (as.numeric(difftime(rightmost_points$datetime[2], rightmost_points$datetime[1], units = "secs")))

  # Rule (3): Rightmost slope cannot be < 1/2 of the steepest ascending slope
  while(!check_last(profile_data)){
    last_ascending_index <- max(which(profile_data$ascending == 1))
 # print(profile_data$slope[last_ascending_index-1])
  #print(max_slope)
  #if (profile_data$slope[last_ascending_index-1] < 0.5 * max_slope){
  #if (rightmost_slope < 0.5 * max_slope) {
   # return(FALSE)
  #print("cutttt")
  # Identify the last ascending index
  #last_ascending_index <- max(which(profile_data$ascending == 1))

  # Set the last point in the ascending segment to not ascending
  profile_data <- profile_data %>%
    dplyr::mutate(
      ascending = dplyr::if_else(
        dplyr::row_number() == last_ascending_index,
        0,
        .data$ascending
      )
    )
  }
  #print("after rule 3")
  #print(profile_data, n = 28)
  #print("post-truncation prof")
  #print(profile_data, n = 28)
  plll<-parallelogram_fit(profile_data)
  #print(plll)
  return(list(profile = profile_data, plll = plll))
}


# Helper function to check if rules are satisfied
check_last <- function(profile_data){
  #print(profile_data, n = 28)
  max_slope <- max(profile_data$slope[profile_data$ascending == 1], na.rm = TRUE)

  # Extract the two rightmost points in the ascending segment
  # Identify the last ascending index
  last_ascending_index <- max(which(profile_data$ascending == 1))
  # print("last ascending")
  # print(last_ascending_index)
  # print(profile_data$slope[last_ascending_index])
  # print("max slope")
  # print(max_slope)
  if (profile_data$slope[last_ascending_index] < 0.5 * max_slope){
    return(FALSE)
  }
  return(TRUE)
}

check_rules <- function(profile_data, max_slope) {
  # Get the ascending segment
  ascending_data <- profile_data %>% dplyr::filter(.data$ascending == 1)

  if (nrow(ascending_data) < 2) {
    return(TRUE)  # If no ascending segment or single point, rules are satisfied
  }

  # Extract the two rightmost points in the ascending segment
  rightmost_points <- ascending_data %>%
    dplyr::arrange(desc(dplyr::row_number())) %>%
    head(2)

  # Calculate the slope of the rightmost segment
  rightmost_slope <- (rightmost_points$melatonin[2] - rightmost_points$melatonin[1]) /
    (as.numeric(difftime(rightmost_points$datetime[2], rightmost_points$datetime[1], units = "secs")))

  # Rule (1a) and (1b): Rightmost slope cannot be zero or negative
  if (rightmost_slope <= 0) {
    return(FALSE)
  }

  # Rule (2a) and (2b): Parallelogram diagonal slope rules
  # TODO commented out 13.12.2024
  # diag1_slope <- (rightmost_points$melatonin[2] - ascending_data$melatonin[1]) /
  #   (as.numeric(difftime(rightmost_points$datetime[2], ascending_data$datetime[1], units = "secs")))
  # diag2_slope <- (ascending_data$melatonin[nrow(ascending_data)] - ascending_data$melatonin[1]) /
  #   (as.numeric(difftime(ascending_data$datetime[nrow(ascending_data)], ascending_data$datetime[1], units = "secs")))
  #
  # diag_ratio <- abs(diag1_slope / diag2_slope)
  # if (diag_ratio < 0.5 || diag_ratio < 0) {
  #   return(FALSE)
  plll<- parallelogram_fit(profile_data)
  #print(plll)
  if(plll$flag){
    return(FALSE)
  }
  return(TRUE)
}

