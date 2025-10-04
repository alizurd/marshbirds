install.packages("spOccupancy")
library(spOccupancy)
library(dplyr)
library(tidyr)
library(janitor)

detection <- read.csv("~/Desktop/marshbirdsoutput/detected.csv", header = TRUE)
detected_clean <- clean_names(detection)

