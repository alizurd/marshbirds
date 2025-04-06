install.packages("caret")
library(caret)

setwd("~/Desktop")
sharp2014 <- read.csv("SHARP_surveyData_2014.csv", stringsAsFactors = FALSE)
View(sharp2014)

colnames(sharp2014) <- tolower(colnames(sharp2014))

sharp2014$tempf <- scale(sharp2014$tempf, center = TRUE, scale = TRUE)
sharp2014$windsp <- scale(sharp2014$windsp, center = TRUE, scale = TRUE)
sharp2014$sky <- scale(sharp2014$sky, center = TRUE, scale = TRUE)
sharp2014$noise <- scale(sharp2014$noise, center = TRUE, scale = TRUE)

#sharp2014$totalcountn <- as.numeric(as.character(sharp2014$totalcountn))

#encoded_tide <- dummyVars("~tide", data = sharp2014)
#data_onehot <- as.data.frame(predict(encoded_tide, newdata = sharp2014))
#sharp2014 <- cbind(encoded_tide, data_onehot)

#sharp2014$totalcountn <- scale(log(sharp2014$totalcountn + 1))

sharpdata <- print(sharp2014)
View(sharpdata)
