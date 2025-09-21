library(unmarked)
library(dplyr)
library(tidyr)
library(janitor)

detection <- read.csv("~/Desktop/marshbirdsoutput/detected.csv", header = TRUE)
detected_clean <- clean_names(detection)

# Filter one species
alphacode <- "CLRA"
df <- detected_clean %>% filter(alphacode == "CLRA")

# create detection matrix
detection_matrix <- df %>%
  select(point_id, year, visit, detection) %>%
  pivot_wider(names_from = visit, values_from = detection, values_fill = 0) %>%
  arrange(point_id, year)

site_covs <- df %>%
  group_by(point_id, year) %>%
  summarize(across(starts_with("hab"), ~first(.)), .groups = "drop")

# unmarked frame
umf <- unmarkedFrameOccu(y = as.matrix(detection_matrix[,3:ncol(detection_matrix)]),
                         siteCovs = site_covs[, -c(1:2)])

# Fit occupancy model
fit <- occu(~1 ~ hab_hm + hab_lm + hab_rd + hab_up + hab_ph, data = umf)
summary(fit)
