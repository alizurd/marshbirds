library(MASS) 
library(janitor)
library(dplyr)
library(lubridate)

# ---------------------------------------------
# Ingest data for GLM
# ---------------------------------------------

sharp_files <- read.csv("~/Desktop/NY_SHARP Survey Data Master File 2011-2024 copy.csv")
sharp_files <- clean_names(sharp_files)

hab_data <- read.csv("~/Desktop/percent_cover_by_plot_jan.csv")

## checking the habitat data for multicollinearity
# hab_data %>%
#   dplyr::select(hm, lm, rd, up, ph) %>%
#   cor(use = "complete.obs") %>%
#   round(2)

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
  group_by(point_id, visit_num, survey_date, time, year, alpha_code) %>%
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
ny_bird_data <- date_time_fix %>%
  filter(
    # state == "NY",
    year %in% c(
                2012,
                2013,
                2014,
                2015,
                2016,
      
                2018,
                2019,
                2021, 
                2022, 
                2023, 
                2024
                ),
    alpha_code %in% c("SALS", "SESP", "CLRA")) # potentially include VIRA and LEBI in future models

# ---------------------------------------------
# Prepare data for GLM
# ---------------------------------------------

# summarise to one row per site per species 
glm_data <- ny_bird_data %>%
  group_by(point_id, alpha_code) %>%
  summarise(
    mean_count = mean(count, na.rm = TRUE),
    total_count = sum(count, na.rm = TRUE),
    n_visits = n(),
    .groups = "drop"
  ) %>%
  left_join(hab_data, by = "point_id") %>%
  filter(!is.na(hm))  # drop sites without habitat data

# # scale habitat covariates
# glm_data <- glm_data %>%
#   mutate(across(c(hm, lm, rd, up, ph), scale))

# ---------------------------------------------
# Run a separate model for each species
# ---------------------------------------------

for (sp in c("CLRA", "SALS", "SESP")) {
  
  cat("\n========================================\n")
  cat("Species:", sp, "\n")
  cat("========================================\n")
  
  sp_data <- glm_data %>% filter(alpha_code == sp)
  
  model <- glm.nb(
    total_count ~ hm + lm + up + ph + offset(log(n_visits)),
    data = sp_data
  )
  
  cat("\n--- Summary ---\n")
  print(summary(model))
  
  assign(paste0("glm_", sp), model)
}
