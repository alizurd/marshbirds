library(janitor)
library(dplyr)

setwd("~/Desktop/marshbirdsoutput/round_8")
richness <- read.csv("richness.csv")
richness <- clean_names(richness)

# --------------------
#   richness model
# --------------------

# glm model with vegetation predictors
model <- glm(richness ~ hm + lm + ow + ph + up + rd, 
             data = richness, 
             family = "poisson")

summary(model)

# ----------------------------------------------

msi <- read.csv("msi.csv")
msi <- clean_names(msi)
View(msi)

# --------------------
#   msi model
# --------------------

# glm model with vegetation predictors

# msi$avg_msi <- as.numeric(msi$avg_msi)

model <- glm(avg_msi ~ hm + lm + ow + ph + up + rd, 
             data = msi, 
             family = "gaussian")

summary(model)

