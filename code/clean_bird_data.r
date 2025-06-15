install.packages("vegan")
install.packages("janitor")
library(janitor)
library(vegan)
library(ggplot2)

setwd("/Users/alyssabueno/Desktop/sharp")
sharp <- read.csv("SHARP_surveyData_2014.csv", stringsAsFactors = FALSE)
View(sharp)

# remove NAs
sharp[is.na(sharp)] <- 0

# remove spaces and apply snakecase to col names
sharp <- clean_names(sharp)

str(sharp) # check column types

# make columns numeric
sharp$total_count_n <- as.numeric(sharp$total_count_n)
sharp$tempf <- as.numeric(sharp$tempf)
sharp$windsp <- as.numeric(sharp$windsp)
sharp$sky <- as.numeric(sharp$sky)
sharp$noise <- as.numeric(sharp$noise)
sharp$survey_date <- as.Date(sharp$surveydate)

keep_columns <- c("totalcountn"
                  , "windsp"
                  , "noise"
                  , "distance")

View(keep_columns)

# checks
nrow(sharp)
unique(sharp$surveydate)
unique(sharp$pointid)
sharp <- aggregate(.~alphacode, data=sharp[,keep_columns], FUN = sum)

decostand(sharp, method = "standardize", margin = 2, na.rm = FALSE)

