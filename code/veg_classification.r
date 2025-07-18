library(terra)
library(randomForest)
library(dplyr)
library(caret)
library(sf)
library(tidyr)
library(ggplot2)

# ------------------------------------
# setting up the training data
# ------------------------------------

setwd ("~/Desktop/marshbirdsoutput/round_5")

# load raster and training polygons
r <- rast("training_raster_round_5.tif") # this has the raster data    
training_polygons <- vect("training_polygons_round_5.shp") # this has the classes

# calculate NDVI
ndvi <- (r[[4]] - r[[1]]) / (r[[4]] + r[[1]])
names(ndvi) <- "ndvi"

# calculate brightness
brightness <- (r[[1]] + r[[2]] + r[[3]] + r[[4]]) / 4
names(brightness) <- "brightness"

# calculate ndwi
ndwi <- (r[[4]] - r[[2]]) / (r[[4]] + r[[2]])
names(ndwi) <- "ndwi"

# now create the PCA 

all_for_pca <- c(r, ndvi)  # add ndvi to the raster
vals <- values(all_for_pca)
vals <- vals[complete.cases(vals), ]
pca_pixels <- prcomp(vals, center = TRUE, scale. = TRUE)
pca <- predict(all_for_pca, pca_pixels, index = 1:5)  # for 5 PCs
names(pca) <- paste0("PCA", 1:5)

# stack em up
r_stack <- c(r, ndvi, brightness, ndwi, pca)

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

# check the sampling balance in the classes
# print("Class distribution in training data:")
# print(table(training_data$class))

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

# -----------------------------
# classifying the plots
# -----------------------------



# load the new, unlabeled raster
prediction_raster <- rast("~/Desktop/marshbirdsoutput/round_4/prediction_raster_round_4.tif")

# calculate NDVI for prediction raster
ndvi_pred <- (prediction_raster[[4]] - prediction_raster[[1]]) / 
  (prediction_raster[[4]] + prediction_raster[[1]])
names(ndvi_pred) <- "ndvi"

# calculation NDWI
ndwi_pred <- (prediction_raster[[4]] - prediction_raster[[2]]) / 
  (prediction_raster[[4]] + prediction_raster[[2]])
names(ndwi_pred) <- "ndwi"

# Calculate brightness
brightness_pred <- (prediction_raster[[1]] + prediction_raster[[2]] + 
                      prediction_raster[[3]] + prediction_raster[[4]]) / 4
names(brightness_pred) <- "brightness"

# Apply the pca
all_for_pca_pred <- c(prediction_raster, ndvi_pred)
names(all_for_pca_pred) <- names(all_for_pca)  # ensure exact match
pca_pred <- predict(all_for_pca_pred, pca_pixels, index = 1:5)
names(pca_pred) <- paste0("PCA", 1:5)

# stack it
prediction_stack <- c(prediction_raster, ndvi_pred, ndwi_pred, brightness_pred, pca_pred)

# rename to match training names
names(prediction_stack) <- names(r_stack)

# apply the model to predict classes
classified <- predict(prediction_stack, rf_model, na.rm = TRUE)

# Define custom colors
colors <- c(
  "#a6d96a",  # hm  → light green
  "#1a9641",  # lm  → dark green
  "#8c510a",  # md  → brown
  "#3288bd",  # ow  → blue
  "#f781bf",  # ph  → pink
  "#636363",  # rd  → gray
  "#762a83"   # up  → purple
)


# Plot with colors
plot(classified, col = colors)

# save the output tif
writeRaster(classified, "~/Desktop/marshbirdsoutput/round_5/classified_output.tif", overwrite = TRUE)

rf_model$importance


# -------------------------------------------------------------------
# feature class importance and visualizations
# -------------------------------------------------------------------


ggplot(training_data, aes(x = PCA3, y = ndvi, color = class)) +
  geom_point(alpha = 0.3) +
  theme_minimal()

ggplot(training_data, aes(x = PCA1, fill = class)) +
  geom_density(alpha = 0.4) +
  theme_minimal()

imp <- as.data.frame(rf_model$importance)
imp$feature <- rownames(imp)

ggplot(imp, aes(x = reorder(feature, MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Random Forest Variable Importance",
       x = "Feature", y = "Mean Decrease Gini") +
  theme_minimal()




# -------------------------------------------------------------------
# calculating class percentages in each plot and exporting as csv
# -------------------------------------------------------------------

# load shp first
plots <- st_read("200m_buffer_clip.shp") 
View(plots)

# convert to terra format for extract()
plots <- vect(plots)

# extract raster values per plot
values <- terra::extract(classified, plots, df = TRUE)

# calculate percent cover per plot
percents <- function(file) {
  r <- values(file)
  freq_table <- as.data.frame(freq(r))  # returns 'value' and 'count'
  
  total_pixels <- sum(freq_table$count, na.rm = TRUE)
  freq_table$percent <- freq_table$count / total_pixels
  
  freq_table$plot_id <- tools::file_path_sans_ext(basename(file))
  return(freq_table)
}

# recreate percents as a dataframe to join
percents <- values %>%
  group_by(ID, class) %>%
  summarize(pixel_count = n(), .groups = "drop") %>%
  group_by(ID) %>%
  mutate(percent = pixel_count / sum(pixel_count)) %>%
  ungroup()

class_labels <- data.frame(
  layer = c(1, 2, 3, 4, 5, 6),  # adjust to match your class codes
  class = c("hm", "lm", "ow", "ph", "rd", "up")
)

class_labels <- dplyr::left_join(percents, class_labels, by = "class")

# pivot to wide format
percent_cover <- class_labels %>%
  select(ID, class, percent) %>%
  pivot_wider(names_from = class, values_from = percent, values_fill = 0) %>%
  rename(plot_id = ID) %>%
  mutate(across(where(is.numeric), \(x) round(x, 2)))

# save to CSV
write.csv(percent_cover, "percent_cover_by_plot.csv", row.names = FALSE)

View(percent_cover)
