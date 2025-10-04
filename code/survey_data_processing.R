library(dplyr)
library(janitor)
library(spOccupancy)

# import data
habitat_perc_cover <- read.csv("~/Desktop/sharp_data/percent_cover_by_plot.csv")

years <- 2011:2014
sharp_files <- sprintf("~/Desktop/sharp_data/working_data/SHARP_surveyData_%d.csv", years)

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

# change the time from hh:mm or h:mm to minutes
full_data <- full_data_raw %>%
  mutate(minutes = sapply(survey_time, function(x) {
    parts <- strsplit(x, ":")[[1]]
    as.numeric(parts[1]) * 60 + as.numeric(parts[2])
  }))

#filter out NY state and 3 species of interest
ny_bird_data <- full_data %>%
  filter(state == "NY",
         alpha_code %in% c("SALS", "SESP", "CLRA"))

hab_data <- habitat_perc_cover %>%
  select(point_id, hm, lm, md, ow, ph, rd, sd, ph, up)

# Join survey data to bird data and filter out null hb rows (indicating we did not have data for those rows)
# based off 2022 survey - potentially look into adding them back in, but need to confirm with liz whether they can provide recent data, which should be prioritized
filtered_data_all <- ny_bird_data %>%
  left_join(hab_data, by = "point_id") %>%
  filter(!is.na(hm)) %>%   # keep only rows with habitat data
  filter(!is.na(lm),      
         !is.na(rd),
         !is.na(up),
         !is.na(ph)) %>%
  filter(!is.na(temp_f),    # filter rows with missing temp and survey time
         !is.na(survey_time)) %>%
  filter(visit_num %in% c("1", "2"))

# convert total_count_n to presence/ absence
filtered_data <- filtered_data_all %>%
  mutate(presence = as.integer(total_count_n > 0)) %>%
  select(-total_count_n)

# View(filtered_data)

# group and summarize to remove any dupes
species_data <- filtered_data %>%
  group_by(point_id, visit_num, survey_date, alpha_code) %>%
  summarize(
    presence = max(presence, na.rm = TRUE),
    # if there are multiple covariate values, take the first one
    hm = first(hm),
    lm = first(lm),
    rd = first(rd),
    up = first(up),
    ph = first(ph),
    observer = first(observer),
    tide = first(tide),
    survey_time = first(survey_time),
    temp_f = first(temp_f),
    sky = first(sky),
    wind_sp = first(wind_sp),
    noise = first(noise),
    .groups = "drop"
  )

#pivot into wide data format
wide_data <- species_data %>%
  select(point_id, visit_num, survey_date, alpha_code, presence) %>%
  pivot_wider(
    names_from = alpha_code, 
    values_from = presence, 
    values_fill = 0 # that might be the issue of missing data
  )

View(wide_data)
# group site level covariates
site_covs <- species_data %>%
  select(point_id, hm, lm, rd, up, ph) %>%
  distinct()

# group survey level covariates
det_covs <- species_data %>%
  select(point_id, visit_num, survey_date, observer, tide, 
         survey_time, temp_f, sky, wind_sp, noise) %>%
  distinct()

# -----------------------------
# create the array 
# -----------------------------

# Create detection array function
create_detection_array <- function(data) {
  species <- c("SALS", "SESP", "CLRA")
  sites <- sort(unique(data$point_id))
  visits <- sort(unique(data$survey_date))
  
  # Initialize array: species x sites x visits
  y_array <- array(
    NA, 
    dim = c(length(species), length(sites), length(visits)),
    dimnames = list(species = species, site = sites, visit = visits)
  )
  
  # Fill array
  for (i in seq_len(nrow(data))) {
    species_idx <- match(data$alpha_code[i], species)
    site_idx <- match(data$point_id[i], sites)
    visit_idx <- match(data$survey_date[i], visits)
    y_array[species_idx, site_idx, visit_idx] <- data$presence[i]
  }
  
  return(y_array)
}

# Create detection covariate matrix (sites x visits)
# Create detection covariate matrix (sites x visits)
create_det_matrix <- function(data, covariate_name) {
  sites <- sort(unique(species_data$point_id))
  visits <- sort(unique(species_data$survey_date))
  
  det_matrix <- matrix(
    NA, 
    nrow = length(sites), 
    ncol = length(visits),
    dimnames = list(site = sites, visit = visits)
  )
  
  for (i in seq_len(nrow(data))) {
    site_idx <- match(data$point_id[i], sites)
    visit_idx <- match(data$survey_date[i], visits)
    # Convert to numeric if needed
    det_matrix[site_idx, visit_idx] <- as.numeric(data[[covariate_name]][i])
  }
  
  # Fill NAs with overall mean
  overall_mean <- mean(det_matrix, na.rm = TRUE)
  det_matrix[is.na(det_matrix)] <- overall_mean
  
  return(det_matrix)
}

# Create array and covariates
y_array <- create_detection_array(species_data)
det_covs_list <- list(
  temp = create_det_matrix(det_covs, "temp_f"),
  time = create_det_matrix(det_covs, "survey_time")
)

# Remove sites with all NA detection histories
sites_with_data <- apply(y_array, 2, function(site_data) {
  any(!is.na(site_data))  # Check if any species/visit has data
})

cat("Sites with data:", sum(sites_with_data), "out of", length(sites_with_data), "\n")

# Filter everything to sites with data
y_array <- y_array[, sites_with_data, ]
det_covs_list$temp <- det_covs_list$temp[sites_with_data, ]
det_covs_list$time <- det_covs_list$time[sites_with_data, ]

sites_to_keep <- dimnames(y_array)$site
site_covs_filtered <- site_covs %>%
  filter(point_id %in% sites_to_keep)

# Build data list
data_list <- list(
  y = y_array,
  occ.covs = site_covs_filtered %>% select(-point_id),
  det.covs = det_covs_list
)

# Verify
cat("\nFinal data check:\n")
cat("y dimensions:", dim(y_array), "\n")
cat("NAs in temp:", sum(is.na(det_covs_list$temp)), "\n")
cat("NAs in time:", sum(is.na(det_covs_list$time)), "\n")
cat("Sites in y:", dim(y_array)[2], "\n")
cat("Sites in occ.covs:", nrow(data_list$occ.covs), "\n")

# Run model
out <- msPGOcc(
  occ.formula = ~ hm + lm,
  det.formula = ~ temp + time,
  data = data_list,
  n.samples = 10000,
  n.burn = 5000,
  n.thin = 5,
  n.chains = 3,
  verbose = TRUE,
  n.report = 1000
)

# Results
summary(out)
plot(out$beta.comm.samples)
save(out, file = "~/Desktop/marshbirds/code/model_results.RData")