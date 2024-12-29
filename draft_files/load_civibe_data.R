# civibe data preparation
library(readr)
saliva_data <- read_csv("~/Salma/tSCN/Projects/civibe/saliva/saliva_data.csv")
civibe_melatonin<- saliva_data%>%dplyr::select(sample_id, timestamp, melatonin_concentration)
colnames(civibe_melatonin)[which(names(civibe_melatonin) == "timestamp")] <- "datetime"
colnames(civibe_melatonin)[which(names(civibe_melatonin) == "melatonin_concentration")] <- "melatonin"

#load corrected melatonin data (ie CivibeFD212 timeswitched samples)
saliva_data_corrected<- saliva_data
saliva_data_corrected[548,]$timestamp <- as.POSIXlt(saliva_data_corrected[548,]$timestamp-2700) #subtract 45 minutes

civibe_melatonin<- saliva_data_corrected%>%dplyr::select(sample_id, timestamp, melatonin_concentration)
colnames(civibe_melatonin)[which(names(civibe_melatonin) == "timestamp")] <- "datetime"
colnames(civibe_melatonin)[which(names(civibe_melatonin) == "melatonin_concentration")] <- "melatonin"


# Define the IDs you want to filter
ids <- c("200", "201", "202", "204", "205", "206", "207", "209", "210", "211", "212", "213") # Add all the IDs you need

# Use purrr::map to filter for each ID and store in a list
civibe_melatonin_FD_list <- purrr::map(ids, ~ civibe_melatonin %>%
                                   dplyr::filter(stringr::str_detect(sample_id, paste0("^CiViBe_", .x, "_FD"))))

# Name the list elements with "mel" prefix and ID
names(civibe_melatonin_FD_list) <- paste0("civibe_melatonin_FD", ids)

# Access the filtered data for each ID
# print(civibe_melatonin_FD_list$civibe_melatonin_FD200, n=45)
#filtered_data_list[["mel201"]] # Example for ID 201
# civibe_201_FD_mel<-civibe_melatonin %>% dplyr::filter(stringr::str_detect(sample_id, "^CiViBe_201_FD"))



##### split data into IDs and days
# Define the IDs you want to filter
ids <- c("200", "201", "202", "204", "205", "206", "207", "209", "210", "211", "212", "213") # Add all the IDs you need

# Use purrr::map to filter and split data
split_data_list <- purrr::map(ids, ~ {
  # Filter data
  filtered_data <- civibe_melatonin %>%
    dplyr::filter(stringr::str_detect(sample_id, paste0("^CiViBe_", .x, "_FD")))

  # Split the data into two halves
  list(
    day1 = filtered_data[1:23, ], # Adjust the indices as needed
    day2 = filtered_data[24:46, ] # Adjust the indices as needed
  )
})

# Name the list elements with "civibe_melatonin_FD<ID>_day1/day2"
names(split_data_list) <- paste0("civibe_melatonin_FD", ids)

# Flatten the list to name day1 and day2
# Split the list into day1 and day2, and properly name them
civibe_melatonin_FD_daysplit <- purrr::imap(split_data_list, ~ {
  # Create a named list for day1 and day2
  setNames(
    list(.x$day1, .x$day2),
    c(paste0(.y, "_day1"), paste0(.y, "_day2"))
  )
}) %>%
  purrr::flatten()


######## add decimalhours column to daysplit data

# Apply the transformation to all elements of the list
civibe_melatonin_FD_dh_daysplit <- civibe_melatonin_FD_daysplit %>%
  purrr::map(~ .x %>%
        dplyr::mutate(decimalhours = posixct_to_decimal(datetime, datetime)))

civibe_melatonin_FD_dh_daysplit <- civibe_melatonin_FD_dh_daysplit %>%
  purrr::map(~ .x %>%
        dplyr::mutate(decimalhours = round(decimalhours, 5)))  # Round to 5 decimal points
# Access the split data
#civibe_melatonin_FD_daysplit[["civibe_melatonin_FD201_day1"]] # Example for day1 of ID 201
#civibe_melatonin_FD_daysplit[["civibe_melatonin_FD201_day2"]] # Example for day2 of ID 201

####### plot sample data ########

# plot by day

# Create a list to store plots
raw_data_plots_daysplit <- purrr::imap(civibe_melatonin_FD_daysplit, ~ {
  # Generate the plot using plot_profile
  plot<- plot_profile(.x, show_segments = FALSE, show_parallelogram = FALSE, show_roi = FALSE, show_dlmoIP = FALSE, show_fit = FALSE, show_roi_heatmap = FALSE, show_roi_small = FALSE, show_roi_big = FALSE)
  # Return the plot (the name will be assigned later)
  plot
})

# Name the list elements to match the data
names(raw_data_plots_daysplit) <- names(civibe_melatonin_FD_daysplit)

# Access the plots
# raw_data_plots_daysplit[["civibe_melatonin_FD201_day1"]] # Example for ID 201, Day 1


# plot entire protocol (day1+day2)

raw_data_plots <- purrr::imap(civibe_melatonin_FD_list, ~ {
  # Generate the plot using plot_profile
  plot<- plot_profile(.x, show_segments = FALSE, show_parallelogram = FALSE, show_roi = FALSE, show_dlmoIP = FALSE, show_fit = FALSE, show_roi_heatmap = FALSE, show_roi_small = FALSE, show_roi_big = FALSE)
  # Return the plot (the name will be assigned later)
  plot
})

# Name the list elements to match the data
names(raw_data_plots) <- names(civibe_melatonin_FD_list)

# Access the plots
# raw_data_plots[["civibe_melatonin_FD201"]] # Example for ID 201


########### calculate dlmo ##############
# Define a list to store results
dlmo_results <- purrr::imap(civibe_melatonin_FD_daysplit, ~ {
  # Calculate dlmo for the current data
  dlmo_value <- calculate_dlmo(.x)

  # Extract the ID and day from the list name
  name_parts <- stringr::str_match(.y, "civibe_melatonin_FD(\\d+)_day(\\d+)")
  id <- name_parts[2]
  day <- name_parts[3]

  # Save the output with dynamic naming
  assign(paste0("dlmo_FD_", id, "_day", day), dlmo_value, envir = .GlobalEnv)

  # Optionally return the value in the list for review
  dlmo_value
})



