library(dplyr)
library(janitor)
library(spOccupancy)
library(tidyr)
library(lubridate)

# ---------------------------------------------
# importing data
# ---------------------------------------------

# import species data
years <- 2011:2014
sharp_files <- sprintf("~/Desktop/sharp_data/working_data/SHARP_surveyData_%d.csv", years)

# import habitat data
hab_data <- habitat_perc_cover <- read.csv("~/Desktop/sharp_data/percent_cover_by_plot.csv")

# ---------------------------------------------
# cleaning data
# ---------------------------------------------

# define survey cols for sharp_files
survey_cols <- c("region_num", "state", "hexagon", "point_id", "patch_id", 
                 "visit_num", "survey_window", "survey_date", "year", "month", 
                 "day", "observer", "tide", "survey_time", "temp_f", "sky", 
                 "wind_sp", "noise", "alpha_code", "total_count_n")

# read csvs, clean col names, and bind them all together 
full_data_raw <- lapply(sharp_files, function(file) {
  read.csv(file) %>%
    clean_names() %>%
    select(all_of(survey_cols))
}) %>%
  bind_rows()

# date must be converted from chr to date
full_data_raw$survey_date <- ymd(full_data_raw$survey_date)
# class(full_data_raw$survey_date)

# time must be converted to integer
# full_data_raw_time <- full_data_raw %>% # first convert the data into hh:mm, then continue onto integer
#   separate(survey_time, into = c("hour", "minute"), sep = ":", convert = TRUE) %>%
#   mutate(survey_time = sprintf("%02d:%02d", hour, minute)) %>%
#   select(-hour, -minute)
# 
# class(full_data_raw$survey_time)
# full_data_raw_time_new$survey_time <- as.integer(full_data_raw_time$survey_time)

# convert time into minutes since midnight and convert date into julian date
# date_time_fix <- full_data_raw %>%
#   mutate(
#     date = yday(survey_date), # convert to julian
#     time = period_to_seconds(hm(survey_time)) / 60, # minutes since midnight
#     total_count_n = as.numeric(total_count_n)
#   ) %>%
#   group_by(point_id, visit_num, date, time, alpha_code, year) %>%
#   summarise(
#     across(
#       c(patch_id, survey_window, observer, tide, temp_f, wind_sp, sky, noise, state),
#       ~ first(.x), # take the first value within group
#       .names = "{.col}"
#     ),
#     count = sum(total_count_n, na.rm = TRUE), # sum the count of individuals
#     .groups = "drop"
#   )

# just converting just the time here, not the date
date_time_fix <- full_data_raw %>%
  mutate(
    time = period_to_seconds(hm(survey_time)) / 60, # minutes since midnight
    total_count_n = as.numeric(total_count_n),
    presence = as.integer(total_count_n > 0)
  ) %>%
  group_by(point_id, visit_num, survey_date, time, alpha_code) %>%
  summarise(
    across(
      c(patch_id, survey_window, observer, tide, temp_f, wind_sp, sky, noise, state, presence),
      ~ first(.x), # take the first value within group
      .names = "{.col}"
    ),
    count = sum(total_count_n, na.rm = TRUE), # sum the count of individuals
    .groups = "drop"
  )

# check to see what's causing the na coercion
# full_data_raw %>%
#   filter(!grepl("^[0-9.]+$", total_count_n)) %>%
#   distinct(total_count_n)

# ---------------------------------------------
# filtering data to state and species of interest
# ---------------------------------------------

#filter out NY state and 3 species of interest
ny_bird_data <- date_time_fix %>%
  filter(state == "NY",
         alpha_code %in% c("SALS", "SESP", "CLRA")) # potentially include VIRA and LEBI in future models

# ---------------------------------------------
# join habitat data with survey data
# ---------------------------------------------

full_data <- ny_bird_data %>%
  left_join(hab_data, by = "point_id") %>%
  filter(!is.na(hm)) 

# ----------------------------------------------------------------------
# Create detection array function
# ----------------------------------------------------------------------
create_detection_array <- function(data) {
  species <- sort(unique(data$alpha_code))  # Use 'data' not 'full_data'
  sites <- sort(unique(data$point_id))
  dates <- sort(unique(data$survey_date))
  
  y_array <- array(
    NA,
    dim = c(length(species), length(sites), length(dates)),
    dimnames = list(species = species, site = sites, date = dates)
  )
  
  for (i in seq_len(nrow(data))) {
    species_idx <- match(data$alpha_code[i], species)
    site_idx <- match(data$point_id[i], sites)
    date_idx <- match(data$survey_date[i], dates)  # FIXED: was site_idx
    
    y_array[species_idx, site_idx, date_idx] <- data$presence[i]
  }
  
  return(y_array)
}

# ----------------------------------------------------------------------
# Create detection covariate matrix function
# ----------------------------------------------------------------------
create_det_matrix <- function(data, covariate_name) {
  sites <- sort(unique(data$point_id))  # Use 'data' not 'full_data'
  dates <- sort(unique(data$survey_date))
  
  det_matrix <- matrix(
    NA,
    nrow = length(sites),
    ncol = length(dates),
    dimnames = list(site = sites, date = dates)
  )
  
  for (i in seq_len(nrow(data))) {
    site_idx <- match(data$point_id[i], sites)
    date_idx <- match(data$survey_date[i], dates)
    
    value <- data[[covariate_name]][i]
    
    if (is.character(value) || is.factor(value)) {
      value <- as.numeric(as.factor(value))
    } else {
      value <- as.numeric(value)
    }
    
    det_matrix[site_idx, date_idx] <- value
  }  # FIXED: closing brace here
  
  # FIXED: These lines outside the loop
  overall_mean <- mean(det_matrix, na.rm = TRUE)
  det_matrix[is.na(det_matrix)] <- overall_mean
  
  return(det_matrix)
}

# ----------------------------------------------------------------------
# Create arrays and matrices - DO THIS FIRST
# ----------------------------------------------------------------------
y_array <- create_detection_array(full_data)

# Site covariates
occ_covs <- full_data %>%
  select(point_id, hm, lm, rd, up, ph) %>%
  distinct() %>%
  arrange(point_id) %>%
  select(-point_id)

# Detection covariates
det_covs_list <- list(
  temp = create_det_matrix(full_data, "temp_f"),
  time = create_det_matrix(full_data, "time"),
  wind = create_det_matrix(full_data, "wind_sp")
)

# ----------------------------------------------------------------------
# Filter to sites with data
# ----------------------------------------------------------------------
sites_with_data <- apply(y_array, 2, function(site_data) {
  any(!is.na(site_data))
})

cat("Sites with data:", sum(sites_with_data), "out of", length(sites_with_data), "\n")

y_array_filtered <- y_array[, sites_with_data, ]
sites_kept <- dimnames(y_array_filtered)$site

occ_covs_filtered <- occ_covs[sites_with_data, ]

det_covs_filtered <- list(
  temp = det_covs_list$temp[sites_with_data, ],
  time = det_covs_list$time[sites_with_data, ],
  wind = det_covs_list$wind[sites_with_data, ]
)

# ----------------------------------------------------------------------
# Run model
# ----------------------------------------------------------------------
data.msom <- list(
  y = y_array_filtered, 
  occ.covs = occ_covs_filtered, 
  det.covs = det_covs_filtered
)

out.msom <- msPGOcc(
  occ.formula = ~ hm + lm + rd + up + ph,
  det.formula = ~ temp + time + wind,
  data = data.msom,
  n.samples = 5000,
  n.burn = 2500,
  n.thin = 2,
  n.chains = 3,
  verbose = TRUE,
  n.report = 500
)

summary(out.msom)


# # ---------------------------------------------
# # convert data structure
# # ---------------------------------------------
# 
# # convert data into long format
# y_long <- full_data %>%
#   group_by (point_id, date, time, visit_num, year) %>%
#   summarise(
#     across(
#       c(patch_id, survey_window, observer, tide, temp_f, wind_sp, sky, noise, state, alpha_code),
#       ~ first(.x), # take the first value within group
#       .names = "{.col}"
#     ),
#     count = sum(count, na.rm = TRUE), # sum the count of individuals
#     .groups = "drop"
#   ) %>%
#   ungroup() %>%
#   glimpse()
# 
# # set what each dimension of the array should contain
# species <- sort(unique(y_long$alpha_code))
# point <- sort(unique(y_long$point_id))
# N <- length(species)
# J <- length(unique(y_long$point_id))
# K <- 2
# y <- array(NA, dim = c(N, J, K))
# 
# # name the dimernstions of y
# dimnames(y)[[1]] <- species
# dimnames(y)[[2]] <- point
# 
# # str(y)
# 
# # create the array with point, date, and visit data
# for (j in 1:J) { # loop through sites
#   for (k in 1:K) { # loop through visits at each site
#     # extract data for each site/visit combo
#     new_df <- y_long %>%
#       filter (point_id == point[j], visit_num == k)
#     # check if there's more than one date for a given combo
#     if (n_distinct(new_df$date > 1)) {
#       # if there is more than one date, just take the data from the first date
#       new_dates <- unique(sort(new_df$date))
#       new_df <- new_df %>%
#         filter(date == new_dates[1])
#     }
#     # if plot j was sampled during visit k, 
#     # new_df should have at least one row of data (basically, one species will be observed)
#     # if not, assume that there was no sample for that visit
#     if (nrow(new_df)>0) {
#       new_species <- which(species %in% new_df$alpha_code)
#       # set value of 1 to species that were observed
#       y[new_species, j, k] <- 1
#       # set value of 0 for all other species
#       y[-new_species, j, k] <- 0
#     }
#   }
# }
# 
# # str(y)
# # apply(y, 1, sum, na.rm = TRUE) #check count per species
# 
# # format detection covariates
# # new_data <- y_long %>%
# #   group_by (point_id, visit_num) %>%
# #   summarize(date = first(unique(survey_date)),
# #             time = first(unique(survey_time)),
# #             observer = first(unique(observer)),
# #             tide = first(unique(tide)),
# #             temp = first(unique(temp_f)),
# #             sky = first(unique(sky)),
# #             wind = first(unique(wind_sp)),
# #             noise = first(unique(noise))) %>%
# #   ungroup() %>%
# #   glimpse()
#   
# # # want to extract julian date
# # # want to extract the time of day as the number of minutes since midnight
# # date <- matrix(NA, nrow = J, ncol = K)
# # time <- matrix(NA, nrow = J, ncol = K)
# # for (j in 1:J) {
# #   for (k in 1:K) { 
# #     new_vals <- new_date %>%
# #       filter(point_id == point[j], visit_num == k) %>%
# #       mutate(date = yday(date),
# #              time = period_to_seconds(hm(time)/60)) %>%
# #       select(date, time) %>%
# #       arrange(date)
# #     if (nrow(new_vals)>0) {
# #       date[j,k] <- new_vals$date[1]
# #       time[j,k] <- new_vals$time[1]
# #     }
# #       }
# # }
# # 
# # glimpse(new_vals)
# 
# # str(date)
# # str(time)
# 
# # -----------------------------------------
# # adding in the covariates
# # -----------------------------------------
# 
# # points need to be in order from the y array
# points_in_order <- dimnames(new_df)$point
# 
# # group site level covariates
# occ_covs <- hab_data %>%
#   filter(point_id %in% points_in_order) %>%
#   arrange(factor(point_id, levels = points_in_order)) %>%  # Sort to match y array
#   select(hm, lm, rd, up, ph)
# 
# # ------------------------------------------
# # create the matrix for det_covs
# # ------------------------------------------
# 
# # Updated function to handle categorical variables
# create_det_matrix <- function(data, covariate_name) {
#   points <- sort(unique(data$point_id))
#   visits <- sort(unique(data$date))
#   
#   det_matrix <- matrix(
#     NA, 
#     nrow = length(points), 
#     ncol = length(visits),
#     dimnames = list(point = points, visit = visits)
#   )
#   
#   for (i in seq_len(nrow(data))) {
#     site_idx <- match(data$point_id[i], points)
#     visit_idx <- match(data$date[i], visits)
#     
#     value <- data[[covariate_name]][i]
#     
#     # Convert categorical to numeric if needed
#     if (is.character(value) || is.factor(value)) {
#       value <- as.numeric(as.factor(value))  # Convert categories to 1, 2, 3, etc.
#     } else {
#       value <- as.numeric(value)
#     }
#     
#     det_matrix[site_idx, visit_idx] <- value
#   }
#   
#   # Fill NAs with overall mean
#   overall_mean <- mean(det_matrix, na.rm = TRUE)
#   det_matrix[is.na(det_matrix)] <- overall_mean
#   
#   return(det_matrix)
# }
# 
# # NOW create the detection covariates list
# det_covs <- det_covs %>%
#   mutate(
#     tide_numeric = case_when(
#       tide == "Low" ~ 1,
#       tide == "Low/Rising" ~ 2,
#       tide == "High/Rising" ~ 3,
#       tide == "High" ~ 4,
#       tide == "High/Falling" ~ 5,
#       tide == "Low/Falling" ~ 6,
#       TRUE ~ NA_real_
#     )
#   )
# 
# # Then use tide_numeric in your covariates list
# det_covs_list <- list(
#   temp = create_det_matrix(det_covs, "temp_f"),
#   time = create_det_matrix(det_covs, "time"),
#   tide = create_det_matrix(det_covs, "tide_numeric"),  # Use numeric version
#   sky = create_det_matrix(det_covs, "sky"),
#   wind = create_det_matrix(det_covs, "wind_sp"),
#   noise = create_det_matrix(det_covs, "noise")
# )
# 
# # group survey level covariates
# # det_covs <- y_long %>%
# #   select(point_id, visit_num, date, observer, tide, 
# #          time, temp_f, sky, wind_sp, noise) %>%
# #   distinct()
# # 
# # det_covs <- list(
# #   temp = create_det_matrix(det_covs, "temp_f"),
# #   time = create_det_matrix(det_covs, "time"),
# #   tide = create_det_matrix(det_covs, "tide"),
# #   sky = create_det_matrix(det_covs, "sky"),
# #   wind = create_det_matrix(det_covs, "wind_sp"),
# #   noise = create_det_matrix(det_covs, "noise")
# # )
# 
# # Check what new_df looks like
# cat("new_df class:", class(new_df), "\n")
# cat("new_df dimensions:", dim(new_df), "\n")
# str(new_df)
# 
# # If it's a dataframe, you need to convert it to a 3D array
# # Use the create_detection_array function we made earlier
# 
# # Make sure this function is defined:
# create_detection_array <- function(data) {
#   species <- sort(unique(data$alpha_code))
#   sites <- sort(unique(data$point_id))
#   visits <- sort(unique(data$date))
#   
#   y_array <- array(
#     NA, 
#     dim = c(length(species), length(sites), length(visits)),
#     dimnames = list(species = species, site = sites, visit = visits)
#   )
#   
#   for (i in seq_len(nrow(data))) {
#     species_idx <- match(data$alpha_code[i], species)
#     site_idx <- match(data$point_id[i], sites)
#     visit_idx <- match(data$date[i], visits)
#     y_array[species_idx, site_idx, visit_idx] <- data$presence[i]
#   }
#   
#   return(y_array)
# }
# 
# # Create the array from your long-format data
# # Replace with whatever your long data is called (y_long, species_data, date_time_fix filtered)
# y_array <- create_detection_array(y_long)
# 
# # Verify it's correct
# cat("y_array dimensions:", dim(y_array), "\n")
# cat("Is it an array?", is.array(y_array), "\n")
# 
# # Now use y_array in the data list
# data.msom <- list(
#   y = y_array,  # Use y_array, not new_df
#   occ.covs = occ_covs,
#   det.covs = det_covs_list
# )
# 
# # survey covariates aka habitat data
# # data into list object
# data.msom <- list(y = new_df, 
#                   occ.covs = occ_covs, 
#                   det.covs = det_covs)
# 
# # build the model
# out.msom <- msPGOcc(  # lowercase 'm', not 'M'
#   occ.formula = ~ hm + lm + rd + up + ph,  # Formula notation, not the object name
#   det.formula = ~ temp + time + wind,      # Formula notation, not the object name
#   data = data.msom,
#   n.samples = 5000,     # Use n.samples, not n.batch/batch.length
#   n.burn = 2500,
#   n.thin = 2,
#   n.chains = 3,
#   verbose = TRUE,
#   n.report = 500
# )
# 
# # 5. Summary
# summary(out.msom)
# 
