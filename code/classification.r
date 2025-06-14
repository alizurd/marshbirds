library(terra)
library(randomForest)
library(dplyr)

setwd ("~/Desktop/r_inputs")

# load raster and training polygons
r <- rast("training_raster.tif")           # 4-band raster
training_polygons <- vect("training_polygons.shp")  # polygons with class labels

extracted <- extract(r, training_polygons, df = TRUE)

# attach class labels
extracted$class <- training_polygons$class[extracted$ID]

# sample a fixed number per class
training_data <- extracted %>%
  group_by(class) %>%
  sample_n(1000) %>%
  ungroup()

# train model
rf_model <- randomForest(factor(class) ~ ., data = training_data[, -1], ntree = 500)

print(rf_model)

# load the new, unlabeled raster (with same 4 bands as training raster)
prediction_raster <- rast("prediction_raster.tif")

names(prediction_raster)
names(prediction_raster) <- c("training_raster_1", "training_raster_2", "training_raster_3", "training_raster_4")


# apply the model to predict classes
predicted <- predict(prediction_raster, rf_model, na.rm = TRUE)

# save the output classified map
writeRaster(predicted, "classified_output.tif", overwrite = TRUE)

# optional: visualize
plot(predicted)
