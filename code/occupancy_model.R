library(unmarked)
library(dplyr)
library(tidyr)
library(janitor)

detection <- read.csv("~/Desktop/marshbirdsoutput/detected.csv", header = TRUE)
detected_clean <- clean_names(detection)

# Filter one species
alphacode <- "CLRA"
df <- detected_clean %>% filter(alphacode == "CLRA")

# Create detection matrix: sites x visits
detection_matrix <- df %>%
  select(point_id, year, visit, detection) %>%
  pivot_wider(names_from = visit, values_from = detection, values_fill = 0) %>%
  arrange(point_id, year)

# Site covariates
site_covs <- df %>%
  group_by(point_id, year) %>%
  summarize(across(starts_with("habitat"), ~first(.)), .groups = "drop")

# unmarked frame
umf <- unmarkedFrameOccu(y = as.matrix(detection_matrix[,3:ncol(detection_matrix)]),
                         siteCovs = site_covs[, -c(1:2)])

# Fit occupancy model
fit <- occu(~1 ~ habitat1 + habitat2, data = umf)
summary(fit)
