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

# Run single-species models
for (sp in 1:3) {
  species_name <- dimnames(y_array_filtered)$species[sp]
  
  # Get data for this species only
  y_single <- y_array_filtered[sp, , ]
  
  # Remove sites with no observations for THIS species
  sites_observed <- apply(y_single, 1, function(x) any(!is.na(x)))
  
  data_single <- list(
    y = y_single[sites_observed, ],
    occ.covs = occ_covs_filtered[sites_observed, ],
    det.covs = lapply(det_covs_filtered, function(x) x[sites_observed, ])
  )
  
  # Run model
  out <- PGOcc(  # Note: PGOcc for single species, not msPGOcc
    occ.formula = ~ hm + lm + rd + up + ph,
    det.formula = ~ temp + time + wind,
    data = data_single,
    n.samples = 20000,
    n.burn = 10000,
    n.thin = 10,
    n.chains = 3,
    verbose = TRUE
  )
  
  assign(paste0("out_", species_name), out)
}

summary(out)
