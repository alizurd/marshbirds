library(terra)
library(randomForest)
library(dplyr)
library(caret)
library(sf)
library(tidyr)
library(ggplot2)
library(tinytex)
library(xfun)

setwd("/Users/alyssabueno/Desktop/marshbirdsoutput/round_6")

# load raster and training polygons
naip <- rast("training_raster_round_6.tif") # this has the raster data    
training_polygons <- vect("training_polygons_round_6.shp") # this has the classes

names(naip) <- paste0("naip", 1:4) # change name of naip bands layer

# calculate NDVI
ndvi <- (naip[[4]] - naip[[1]]) / (naip[[4]] + naip[[1]])
names(ndvi) <- "ndvi"

# calculate brightness
brightness <- (naip[[1]] + naip[[2]] + naip[[3]] + naip[[4]]) / 4
names(brightness) <- "brightness"

# calculate ndwi
ndwi <- (naip[[4]] - naip[[2]]) / (naip[[4]] + naip[[2]])
names(ndwi) <- "ndwi"

# now create the PCA 

all_for_pca <- c(naip, ndvi)  # add ndvi to the raster
vals <- values(all_for_pca)
vals <- vals[complete.cases(vals), ]
pca_pixels <- prcomp(vals, center = TRUE, scale. = TRUE)
pca <- predict(all_for_pca, pca_pixels, index = 1:5)  # for 5 PCs
names(pca) <- paste0("PCA", 1:5)

# stack em up
r_stack <- c(pca)

# extract raster values for training polygons
extracted <- terra::extract(r_stack, training_polygons, df = TRUE)

# turn class into a factor
training_polygons$class <- as.factor(training_polygons$class)

# Convert training_polygons to dataframe to get the class labels
poly_df <- as.data.frame(training_polygons)
poly_df$ID <- 1:nrow(poly_df)  # Add ID column to match extract output

# Join the class labels
extracted <- extracted %>%
  left_join(poly_df[, c("ID", "class")], by = "ID")

# clean data by removing rows with na values and remove ID column
extracted_clean <- na.omit(extracted[, -1])  # remove ID and NAs

# sample 2000 per class
training_data <- extracted_clean %>%
  group_by(class) %>%
  sample_n(min(2000, n())) %>%  # Use min() to handle small classes
  ungroup()

# turn class column into factor
training_data$class <- factor(training_data$class)

# split training and test data
set.seed(342)
idx <- sample(seq_len(nrow(training_data)), size = 0.8 * nrow(training_data))
train_set <- training_data[idx, ]
test_set  <- training_data[-idx, ]

# train random forest model
rf_model <- randomForest(class ~ ., 
                         data = train_set, 
                         ntree = 500)

# validation metrics
preds <- predict(rf_model, newdata = test_set)
accuracy <- mean(preds == test_set$class)
cat("Accuracy:", round(accuracy, 4), "\n")

print(rf_model)

# confusion matrix
confusionMatrix(preds, test_set$class)

# directory with your 19 subset rasters
tile_dir <- "~/Desktop/marshbirdsoutput/subset_rasters/"
tile_files <- list.files(tile_dir, pattern = "\\.tif$", full.names = TRUE)

out_dir <- "~/Desktop/marshbirdsoutput/classified_rasters/"
dir.create(out_dir, showWarnings = FALSE)

# define color paletter
colors <- rainbow(length(levels(train_set$class)))

# loop through each subset raster
for (tile in tile_files) {
  tryCatch({
    message("Processing: ", tile)
    
    # load subset tile
    prediction_raster <- rast(tile)
    names(prediction_raster) <- paste0("naip", 1:4)
    
    # compute ONLY what you used in training
    ndvi_tile <- (prediction_raster[[4]] - prediction_raster[[1]]) / 
      (prediction_raster[[4]] + prediction_raster[[1]])
    names(ndvi_tile) <- "ndvi"
    
    # Match training exactly - only naip + ndvi
    all_for_pca_pred <- c(prediction_raster, ndvi_tile)
    names(all_for_pca_pred) <- names(all_for_pca)
    
    # PCA transform
    pca_pred <- terra::predict(all_for_pca_pred, pca_pixels, index = 1:5)
    names(pca_pred) <- paste0("PCA", 1:5)
    
    # match training stack
    prediction_stack <- pca_pred
    names(prediction_stack) <- names(r_stack)
    
    # output filename
    outfile <- file.path(
      out_dir,
      paste0(tools::file_path_sans_ext(basename(tile)), "_classified.tif")
    )
    
    # classify and save
    classified <- terra::predict(
      prediction_stack, rf_model,
      na.rm = TRUE,
      filename = outfile,
      overwrite = TRUE
    )
    
    plot(classified, col = colors, main = basename(tile))
    
    # clean up memory
    rm(prediction_raster, ndvi_tile, all_for_pca_pred, pca_pred, prediction_stack, classified)
    gc()
    
  }, error = function(e) {
    message("Error processing ", tile, ": ", e$message)
  })
}

# list all classified tif files
classified_files <- list.files(classified_dir, pattern = "\\.tif$", full.names = TRUE)

# load all classified rasters
classified_rasters <- lapply(classified_files, rast)

# mosaic them together
# do.call with terra::mosaic combines all rasters
mosaic_raster <- do.call(mosaic, c(classified_rasters, fun = "max"))

# save the mosaic
mosaic_outfile <- "~/Desktop/marshbirdsoutput/mosaic_classified.tif"
writeRaster(mosaic_raster, mosaic_outfile, overwrite = TRUE)

# Define custom colors
colors <- c(
  "#a6d96a",  # hm  
  "#1a9641",  # lm  
  "#8c510a",  # md  
  "#3288bd",  # ow  
  "#fdae61",  # ph  
  "#969696",  # rd  
  "#762a83"   # up  
)



# Plot with colors
plot(classified, col = colors)

# save the output tif
writeRaster(classified, "~/Desktop/marshbirdsoutput/round_6/classified_output.tif", overwrite = TRUE)

rf_model$importance

ggplot(training_data, aes(x = PCA3, y = PCA1, color = class)) +
  geom_point(alpha = 0.3) +
  theme_minimal()

ggplot(training_data, aes(x = PCA1, fill = class)) +
  geom_density(alpha = 0.4, color = NA) +
  theme_minimal()

imp <- as.data.frame(rf_model$importance)
imp$feature <- rownames(imp)

ggplot(imp, aes(x = reorder(feature, MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Random Forest Variable Importance",
       x = "Feature", y = "Mean Decrease Gini") +
  theme_minimal()

