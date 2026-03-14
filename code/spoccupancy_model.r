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
sharp_files <- read.csv("~/Desktop/NY_SHARP Survey Data Master File 2011-2024 copy.csv")
sharp_files <- clean_names(sharp_files)
# View(sharp_files)

# import habitat data
hab_data <- read.csv("~/Desktop/percent_cover_by_plot_jan.csv")

# ---------------------------------------------
# cleaning data
# ---------------------------------------------

# date must be converted from chr to date
sharp_files$survey_date <- mdy(sharp_files$survey_date)
# class(full_data_raw$survey_date)

# create year column from survey date 
sharp_files$year <- year(sharp_files$survey_date)  

# define survey cols for sharp_files
survey_cols <- c("region_num", "state", "hexagon", "point_id", "patch_id", 
                 "visit_num", "survey_date", "year",
                 "observer", "tide", "survey_time", "temp_f", "sky", 
                 "wind_sp", "noise", "alpha_code", "total_count")

# converts nulls to NA to silence any warnings when converting time
sharp_files$survey_time <- ifelse(sharp_files$survey_time == "NULL", NA, sharp_files$survey_time)

# just converting just the time here, not the date
date_time_fix <- sharp_files %>%
  mutate(
    time = period_to_seconds(hm(survey_time)) / 60,
    total_count = as.numeric(total_count),
  ) %>%
  group_by(point_id, visit_num, survey_date, time, alpha_code) %>%
  summarise(
    across(
      c(observer, tide, temp_f, wind_sp, sky, noise),
      ~ first(.x),
      .names = "{.col}"
    ),
    count = sum(total_count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(presence = as.integer(count > 0))  

# ---------------------------------------------
# filtering data to state and species of interest
# ---------------------------------------------

#filter out NY state and 3 species of interest -- update: new files that Liz sent already filter state
# ny_bird_data <- date_time_fix %>%
#   filter(
    # state == "NY",
         # alpha_code %in% c("SALS", "SESP", "CLRA"
                           )) # potentially include VIRA and LEBI in future models

# ---------------------------------------------
# join habitat data with survey data
# ---------------------------------------------

full_data <- date_time_fix %>%
  left_join(hab_data, by = "point_id") %>%
  filter(!is.na(hm)) 


# create a complete grid of all site-visit-species combos
all_visits <- full_data %>%
  distinct(point_id, survey_date)

focal_species <- c("CLRA", "SESP", "SALS")

complete_grid <- all_visits %>%
  crossing(alpha_code = focal_species)

# join actual observations onto the complete grid
focal_species_data <- complete_grid %>%
  left_join(
    full_data %>% filter(alpha_code %in% focal_species),
    by = c("point_id", "survey_date", "alpha_code")
  ) %>%
  mutate(presence = ifelse(is.na(presence), 0L, presence))

# # ----------------------------------------------------------------------
# # filter out the sites that don't have any survey data
# # ----------------------------------------------------------------------
# 
# #filter for sites that are not completely NA
# points_with_data <- full_data %>%
#   group_by(point_id) %>%
#   summarise(contains_data = any(!is.na(presence)), .groups = 'drop') %>% # will return a single true if any observation is not NA
#   filter(contains_data) %>%
#   pull(point_id)
# 
# full_data_clean <- full_data %>%
#   filter(point_id %in% points_with_data)

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
y_array <- create_detection_array(focal_species_data)

# Site covariates
occ_covs <- focal_species_data %>%
  distinct(point_id, hm, lm, rd, ph) %>%
  arrange(point_id) %>%
  dplyr::select(-point_id)

# Detection covariates
det_covs <- list(
  temp = create_det_matrix(focal_species_data, "temp_f"),
  time = create_det_matrix(focal_species_data, "time"),
  wind = create_det_matrix(focal_species_data, "wind_sp"),
  noise = create_det_matrix(focal_species_data, "noise")
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
    occ.formula = ~ hm + lm + rd + ph,
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

## checks

# what % of site-visits actually detected each species?
date_time_fix %>%
  filter(alpha_code %in% c("SALS", "SESP", "CLRA")) %>%
  group_by(alpha_code) %>%
  summarise(
    total_visits = n(),
    detections = sum(presence, na.rm = TRUE),
    detection_rate = mean(presence, na.rm = TRUE)
  )
