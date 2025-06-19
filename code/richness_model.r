library(janitor)
library(dplyr)

setwd("~/Desktop/r_inputs")
bird <- read.csv("SHARP_surveyData_2014.csv")
bird <- clean_names(bird)
View(bird)

# --------------------------
#   preparing bird data
# --------------------------

species_richness <- bird %>%
  group_by(point_id) %>%
  summarize(richness = n_distinct(alpha_code), .groups = 'drop')

# Filter for specific point IDs and rename them
bx_birds <- species_richness %>%
  filter(point_id %in% c("PEL1", "PEL2", "PEL3")) %>%
  mutate(point_id = case_when(
    point_id == "PEL1" ~ "1",
    point_id == "PEL2" ~ "2", 
    point_id == "PEL3" ~ "3",
    TRUE ~ point_id
  )) %>%
  rename(plot_id = point_id)

# ------------------------------
#   adding and joining veg data
# ------------------------------

veg <- read.csv("percent_cover_by_plot.csv")
veg <- clean_names(veg)

# Convert plot_id to numeric for joining
bx_birds$plot_id <- as.numeric(bx_birds$plot_id)

# Join the datasets
data <- bx_birds %>%
  left_join(veg, by = "plot_id")

# --------------------
#   model v1
# --------------------

# Create GLM model with vegetation predictors
model <- glm(richness ~ hm + lm + ow + ph + up + rd, 
             data = data, 
             family = "poisson")

summary(model)
