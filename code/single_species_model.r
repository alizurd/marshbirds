library(dplyr)
library(janitor)
library(spOccupancy)
library(tidyr)
library(lubridate)

# ---------------------------------------------
# importing data
# ---------------------------------------------

# import species data
# years <- 2011:2014
# sharp_files <- sprintf("~/Desktop/sharp_data/working_data/SHARP_surveyData_%d.csv", years)
sharp_files <- read.csv("~/Desktop/NY_SHARP Survey Data Master File 2011-2024 - Copy of NY_SHARP Survey Data Master File 2011-2024 (1).csv")
sharp_files <- clean_names(sharp_files)
View(sharp_files)

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
  filter(!grepl("^[0-9.]+$", total_count_n)) %>%
  distinct(total_count_n)

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
# filter out the sites that don't have any survey data
# ----------------------------------------------------------------------

#filter for sites that are not completely NA
points_with_data <- full_data %>%
  group_by(point_id) %>%
  summarise(contains_data = any(!is.na(presence)), .groups = 'drop') %>% # will return a single true if any observation is not NA
  filter(contains_data) %>%
  pull(point_id)

full_data_clean <- full_data %>%
  filter(point_id %in% points_with_data)

# ----------------------------------------------------------------------
# Create detection array function
# ----------------------------------------------------------------------
create_detection_array <- function(data) {
  species <- sort(unique(data$alpha_code))  
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
    date_idx <- match(data$survey_date[i], dates)  
    
    y_array[species_idx, site_idx, date_idx] <- data$presence[i]
  }
  
  return(y_array)
}

# ----------------------------------------------------------------------
# Create detection covariate matrix function
# ----------------------------------------------------------------------
create_det_matrix <- function(data, covariate_name) {
  sites <- sort(unique(data$point_id)) 
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
  }  
  
  overall_mean <- mean(det_matrix, na.rm = TRUE)
  det_matrix[is.na(det_matrix)] <- overall_mean
  
  return(det_matrix)
}

# ----------------------------------------------------------------------
# Create arrays and matrices 
# ----------------------------------------------------------------------
y_array <- create_detection_array(full_data_clean)

# Site covariates
occ_covs <- full_data_clean %>%
  select(point_id, hm, lm, rd, up, ph) %>%
  distinct() %>%
  arrange(point_id) %>%
  select(-point_id)

# Detection covariates
det_covs <- list(
  temp = create_det_matrix(full_data_clean, "temp_f"),
  time = create_det_matrix(full_data_clean, "time"),
  wind = create_det_matrix(full_data_clean, "wind_sp"),
  noise = create_det_matrix(full_data_clean, "noise")
)


# Run single-species models
for (sp in 1:3) {
  species_name <- dimnames(y_array)$species[sp]
  
  cat("\n========================================\n")
  cat("Running model for species:", species_name, "\n")
  cat("========================================\n\n")
  
  # Get data for this species only
  y <- y_array[sp, , ]
  
  # Remove sites with no observations for THIS species
  sites_observed <- apply(y, 1, function(x) any(!is.na(x)))
  
  data_filtered <- list(
    y = y[sites_observed, ],
    occ.covs = occ_covs[sites_observed, ],
    det.covs = lapply(det_covs, function(x) x[sites_observed, ])
  )
  
  
  # Run model
  out <- PGOcc(  
    occ.formula = ~ hm + lm + rd + up + ph,
    det.formula = ~ temp + noise + wind,
    data = data_filtered,
    n.samples = 20000,
    n.burn = 10000,
    n.thin = 10,
    n.chains = 3,
    verbose = TRUE
  )
  
  assign(paste0("out_", species_name), out)
  
}

# summarize each species 
cat("\n======== CLRA RESULTS ========\n")
summary(out_CLRA)

cat("\n======== SALS RESULTS ========\n")
summary(out_SALS)

cat("\n======== SESP RESULTS ========\n")
summary(out_SESP)
