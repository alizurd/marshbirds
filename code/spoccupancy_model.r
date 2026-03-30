library(dplyr)
library(janitor)
library(spOccupancy)
library(tidyr)
library(lubridate)

# ---------------------------------------------
# importing data
# ---------------------------------------------

sharp_files <- read.csv("~/Desktop/NY_SHARP Survey Data Master File 2011-2024 copy.csv")
sharp_files <- clean_names(sharp_files)

hab_data <- read.csv("~/Desktop/percent_cover_by_plot_jan.csv")

# ---------------------------------------------
# cleaning data
# ---------------------------------------------

sharp_files$survey_date <- mdy(sharp_files$survey_date)
sharp_files$year <- year(sharp_files$survey_date)  
sharp_files$survey_time <- ifelse(sharp_files$survey_time == "NULL", NA, sharp_files$survey_time)

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
# join habitat data with survey data
# ---------------------------------------------

# filter to only sites that have habitat data
sites_with_hab <- hab_data %>% distinct(point_id)

full_data <- date_time_fix %>%
  filter(point_id %in% sites_with_hab$point_id) %>%
  # collapse multiple visits on same date to one row per site-date-species
  # presence = 1 if detected in any visit; covariates take the first visit
  group_by(point_id, survey_date, alpha_code) %>%
  summarise(
    across(c(observer, tide, temp_f, wind_sp, sky, noise, time), ~ first(.x)),
    presence = max(presence, na.rm = TRUE),
    .groups = "drop"
  )

# ---------------------------------------------
# build complete site x date x species grid
# ---------------------------------------------

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
  mutate(presence = ifelse(is.na(presence), 0L, presence)) %>%
  # join hab data here so every row (including zero-fills) gets hab values
  left_join(hab_data, by = "point_id")

# ---------------------------------------------
# filter to sites that have survey data
# ---------------------------------------------

points_with_data <- focal_species_data %>%
  group_by(point_id) %>%
  summarise(contains_data = any(!is.na(presence)), .groups = 'drop') %>%
  filter(contains_data) %>%
  pull(point_id)

full_data_clean <- focal_species_data %>%
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

# site covariates — one row per site
occ_covs <- full_data_clean %>%
  group_by(point_id) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(point_id, hm, lm, ph) %>%
  arrange(point_id)

# verify one row per site before proceeding
stopifnot(nrow(occ_covs) == n_distinct(full_data_clean$point_id))

# detection covariates
det_covs <- list(
  temp  = create_det_matrix(full_data_clean, "temp_f"),
  time  = create_det_matrix(full_data_clean, "time"),
  wind  = create_det_matrix(full_data_clean, "wind_sp"),
  noise = create_det_matrix(full_data_clean, "noise")
)

# ----------------------------------------------------------------------
# Run single-species models
# ----------------------------------------------------------------------

for (sp in 1:3) {
  species_name <- dimnames(y_array)$species[sp]
  
  cat("\n========================================\n")
  cat("Running model for species:", species_name, "\n")
  cat("========================================\n\n")
  
  y <- y_array[sp, , ]
  
  sites_observed <- apply(y, 1, function(x) any(!is.na(x)))
  observed_site_names <- names(sites_observed)[sites_observed]
  
  data_filtered <- list(
    y = y[sites_observed, ],
    occ.covs = occ_covs %>%
      filter(point_id %in% observed_site_names) %>%
      arrange(point_id) %>%
      dplyr::select(-point_id),
    det.covs = lapply(det_covs, function(x) x[sites_observed, ])
  )
  
  out <- PGOcc(  
    occ.formula = ~ hm + lm + ph,
    det.formula = ~ temp 
    + noise 
    # + wind 
    # + time
    ,
    data = data_filtered,
    n.samples = 20000,
    n.burn = 10000,
    n.thin = 10,
    n.chains = 3,
    verbose = TRUE
  )
  
  assign(paste0("out_", species_name), out)
}

# ----------------------------------------------------------------------
# Summarize results
# ----------------------------------------------------------------------

cat("\n======== CLRA RESULTS ========\n")
summary(out_CLRA)

cat("\n======== SALS RESULTS ========\n")
summary(out_SALS)

cat("\n======== SESP RESULTS ========\n")
summary(out_SESP)

# ----------------------------------------------------------------------
# Checks
# ----------------------------------------------------------------------

# detection rates per species (using modeled data only — sites with habitat data)
full_data_clean %>%
  filter(alpha_code %in% focal_species) %>%
  group_by(alpha_code) %>%
  summarise(
    total_visits = n(),
    detections = sum(presence, na.rm = TRUE),
    detection_rate = mean(presence, na.rm = TRUE)
  )
