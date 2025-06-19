library(terra)
library(randomForest)
library(dplyr)
library(caret)
library(sf)
library(tidyr)

# ------------------------------------
# setting up the training data
# ------------------------------------

setwd ("~/Desktop/r_inputs")

# load raster and training polygons
r <- rast("training_raster.tif") # this has the raster data    
training_polygons <- vect("training_polygons.shp") # this has the classes

# calculate NDVI
ndvi <- (r[[4]] - r[[1]]) / (r[[4]] + r[[1]])
names(ndvi) <- "ndvi"

# calculate brightness
brightness <- (r[[1]] + r[[2]] + r[[3]] + r[[4]]) / 4
names(brightness) <- "brightness"

# calculate ndwi
ndwi <- (r[[4]] - r[[2]]) / (r[[4]] + r[[2]])
names(ndwi) <- "ndwi"

# stack em up
r_stack <- c(r, ndvi, brightness, ndwi)

# add r_stack to extracted
extracted <- extract(r_stack, training_polygons, df = TRUE)

# attach class labels
extracted$class <- training_polygons$class[extracted$ID]

# clean data
training_data <- na.omit(extracted[, -1])  # remove ID

# sample a fixed number per class
training_data <- extracted %>%
  group_by(class) %>%
  sample_n(1000) %>%
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

# confusion matrix
confusionMatrix(preds, test_set$class)

# -----------------------------
# classifying the plots
# -----------------------------

# load the new, unlabeled raster
prediction_raster <- rast("prediction_raster.tif")

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

# stack it
prediction_stack <- c(prediction_raster, ndvi_pred, ndwi_pred, brightness_pred)

# rename to match training names but without ID
names(prediction_stack) <- names(r_stack)[names(r_stack) != "ID"]

# apply the model to predict classes
classified <- predict(prediction_stack, rf_model, na.rm = TRUE)

plot(classified)

# save the output tif
writeRaster(classified, "classified_output.tif", overwrite = TRUE)

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

# recreate percents as a dataframe so you can join
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
