install.packages("randomForest")
library(tidyverse)
library(randomForest)

setwd ("/Users/alyssabueno/Desktop")
veg_data <- read.csv("bx_veg_data.csv", stringsAsFactors = FALSE)  
veg_data$class <- as.factor(veg_data$class)

set.seed(123)
train_idx <- sample(nrow(veg_data), 0.8 * nrow(veg_data))
train_data <- veg_data[train_idx, ]
test_data <- veg_data[-train_idx, ]

rf_model <- randomForest(class ~ band_1 + band_2 + band_3, data = train_data, ntree = 500)
print(rf_model)

# Evaluate Performance

predictions <- predict(rf_model, test_data)
confusionMatrix <- table(Predicted = predictions, Actual = test_data$class)
print(confusionMatrix)
