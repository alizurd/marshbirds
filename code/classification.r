library(terra)
library(randomForest)
library(dplyr)
library(caret)

setwd ("~/Desktop/r_inputs")

# load raster and training polygons
r <- rast("training_raster.tif") # this has the raster data    
training_polygons <- vect("training_polygons.shp") # this has the classes

# calculate NDVI
ndvi <- (r[[4]] - r[[1]]) / (r[[4]] + r[[1]])
names(ndvi) <- "ndvi"

# stack em up
r_stack <- c(r, ndvi)

# add r_stack to extracted
extracted <- extract(r_stack, training_polygons, df = TRUE)

# attach class labels
extracted$class <- training_polygons$class[extracted$ID]

# clean data
training_data <- na.omit(extracted[, -1])  # remove ID

# sample a fixed number per class
training_data <- extracted %>%
  group_by(class) %>%
  sample_n(10000) %>%
  ungroup()

# rename the high marsh column
training_data$class[training_data$class == "hm1"] <- "hm"
sum(is.na(training_data))
training_data <- na.omit(training_data)

# turn class column into factor
training_data$class <- factor(training_data$class)

#check the sampling balance in the classes
# table(training_data$class)
# View(training_data)

# split training and test data
set.seed(342)
idx <- sample(seq_len(nrow(training_data)), size = 0.8 * nrow(training_data))
train_set <- training_data[idx, ]
test_set  <- training_data[-idx, ]

# model without ID col
rf_model <- randomForest(class ~ ., 
                         data = train_set[, -1], 
                         ntree = 500)

# validation metrics
preds <- predict(rf_model, newdata = test_set)
mean(preds == test_set$class)

print(rf_model)

# load the new, unlabeled raster
prediction_raster <- rast("prediction_raster.tif")

# calculate NDVI for prediction raster
ndvi_pred <- (prediction_raster[[4]] - prediction_raster[[1]]) / 
  (prediction_raster[[4]] + prediction_raster[[1]])
names(ndvi_pred) <- "ndvi"

# stack it
prediction_stack <- c(prediction_raster, ndvi_pred)

# rename to match training names but without ID
names(prediction_stack) <- names(r_stack)[names(r_stack) != "ID"]

# apply the model to predict classes
predicted <- predict(prediction_stack, rf_model, na.rm = TRUE)

plot(predicted)

# save the output tif
writeRaster(predicted, "classified_output_with_ndvi.tif", overwrite = TRUE)

