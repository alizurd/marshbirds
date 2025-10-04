library(dplyr)
library(janitor)
library(spOccupancy)
library(tidyr)

# ---------------------------------------------
# importing data
# ---------------------------------------------

# import data
habitat_perc_cover <- read.csv("~/Desktop/sharp_data/percent_cover_by_plot.csv")

years <- 2011:2014
sharp_files <- sprintf("~/Desktop/sharp_data/working_data/SHARP_surveyData_%d.csv", years)

# import habitat data
habitat_perc_cover <- read.csv("~/Desktop/sharp_data/percent_cover_by_plot.csv")

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

#filter out NY state and 3 species of interest
ny_bird_data <- full_data_raw %>%
  filter(state == "NY",
         alpha_code %in% c("SALS", "SESP", "CLRA")) # potentially include VIRA and LEBI in future models

# ---------------------------------------------
# convert data structure
# ---------------------------------------------

# convert data into long format
y_long <- ny_bird_data %>%
  group_by (point_id, survey_date, visit_num, alpha_code) %>%
  summarize(total_count_n = n()) %>%
  ungroup() %>%
  glimpse()

# set what each dimension of the array should contain
species <- sort(unique(y_long$alpha_code))
point <- sort(unique(y_long$point_id))
N <- length(species)
J <- length(unique(y_long$point_id))
K <- 2
y <- array(NA, dim = c(N, J, K))

# name the dimernstions of y
dimnames(y)[[1]] <- species
dimnames(y)[[2]] <- point

# str(y)

for (j in 1:J) { # loop through sites
  for (k in 1:K) { # loop through visits at each site
    # extract data for each site/visit combo
    new_df <- y_long %>%
      filter (point_id == point[j], visit_num == k)
    # check if there's more than one date for a given combo
    if (n_distinct(new_df$survey_date > 1)) {
      # if there is more than one date, just take the data from the first datwae
      new_dates <- unique(sort(new_df$survey_date))
      new_df <- new_df %>%
        filter(survey_date == new_dates[1])
    }
    # if plot j was sampled during visit k, 
    # new_df should have at least one row of data (basically, one species will be observed)
    # if not, assume that there was no sample for that visit
    if (nrow(new_df)>0) {
      new_species <- which(species %in% new_df$alpha_code)
      # set value of 1 to species that were observed
      y[new_species, j, k] <- 1
      # set value of 0 for all other species
      y[-new_species, j, k] <- 0
    }
  }
}

# str(y)
# apply(y, 1, sum, na.rm = TRUE) #check count per species

# format detection covariates
new_date <- ny_bird_data %>%
  group_by (point_id, visit_num) %>%
  summarize(date = first(unique(survey_date)),
            time = first(unique(survey_time))) %>%
  ungroup() %>%
  glimpse()

# want to extract julian date
# want to extract the time of day as the number of minutes since midnight
date <- matrix(NA, nrow = J, ncol = K)
time <- matrix(NA, nrow = J, ncol = K)
for (j in 1:J) {
  for (k in 1:K) { 
    new_vals <- new_date %>%
      filter(point_id == point[j], visit_num == k) %>%
      select(date, time) %>%
      arrange(date)
    if (nrow(new_vals)>0) {
      date[j,k] <- new_vals$date[1]
      time[j,k] <- new_vals$time[1]
    }
      }
}

# str(date)
# str(time)

# -----------------------------------------
# adding in the covariates
# -----------------------------------------

# group site level covariates
occ_covs <- hab_data %>%
  select(point_id, hm, lm, rd, up, ph) %>%
  distinct()

# group survey level covariates
# figure out a way to bring in all these other det covariates
det_covs <- species_data %>%
  select(point_id, visit_num, survey_date, observer, tide, 
         survey_time, temp_f, sky, wind_sp, noise) %>%
  distinct()

# detection covariates
det_covs <- list(day = date,
                 time_of_day = time)

str(det_covs)

# survey covariates aka habitat data
# data into list object
data.msom <- list(y = new_df, 
                  occ.covs = occ_covs, 
                  det.covs = det_covs)

# build the model
out.msom <- spMsPGOcc(occ_formula = ~ occ.covs), 
                      det_formula = ~ det.covs, 
                      data = data.msom,
                      n.batch = 10, 
                      batch.length = 25, 
                      cov.model = 'exponential', 
                      NNGP = TRUE, 
                      verbose = FALSE) 
summary(out.msom, level = 'community')

