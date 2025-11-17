library(terra)
library(randomForest)
library(dplyr)
library(caret)
library(sf)
library(tidyr)
library(ggplot2)
library(tinytex)
library(xfun)

setwd("/Users/alyssabueno/Desktop/marshbirdsoutput/round_7")

# load raster and training polygons
naip <- rast("training_raster_round_7.tif") # this has the raster data    
training_polygons <- vect("training_polygons_all.shp") # this has the classes

names(naip) <- paste0("naip", 1:4) # change name of naip bands layer

# calculate NDVI
ndvi <- (naip[[4]] - naip[[1]]) / (naip[[4]] + naip[[1]])
names(ndvi) <- "ndvi"

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

# Join the class labels -- might need to adjust this part since im removing the ID field from the shapefiles??
extracted <- extracted %>%
  left_join(poly_df[, c("ID", "class")], by = "ID")

# clean data by removing rows with na values and remove ID column
extracted_clean <- na.omit(extracted[, -1])  # remove ID and NAs

# sample 4000 per class
training_data <- extracted_clean %>%
  group_by(class) %>%
  sample_n(min(4000, n())) %>%  # Use min() to handle small classes
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

class_colors <- c(
  "hm" = "#a6d96a",  # high marsh
  "lm" = "#1a9641",  # low marsh  
  "md" = "#8c510a",  # mudflat
  "ow" = "#3288bd",  # open water
  "ph" = "#fdae61",  # phragmites
  "rd" = "#969696",  # road
  "up" = "#762a83",   # upland
  "sd" = "#cdae65" # sand
)

# define color paletter
colors <- class_colors[levels(train_set$class)]

# loop through each subset raster
for (tile in tile_files) {
  tryCatch({
    message("Processing: ", tile)
    
    # load subset tile
    prediction_raster <- rast(tile)
    names(prediction_raster) <- paste0("naip", 1:4)
    
    # compute what was used in training
    ndvi_tile <- (prediction_raster[[4]] - prediction_raster[[1]]) / 
      (prediction_raster[[4]] + prediction_raster[[1]])
    names(ndvi_tile) <- "ndvi"
    
    all_for_pca_pred <- c(prediction_raster, ndvi_tile)
    names(all_for_pca_pred) <- names(all_for_pca)
    
    # PCA transform
    pca_pred <- terra::predict(all_for_pca_pred, pca_pixels, index = 1:5)
    names(pca_pred) <- paste0("PCA", 1:5)
    
    # match training stack
    prediction_stack <- pca_pred
    names(prediction_stack) <- names(r_stack)
    
    outfile <- file.path(
      out_dir,
      paste0(tools::file_path_sans_ext(basename(tile)), "_classified.tif")
    )
    
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

# print random forest importance
print("Random Forest Variable Importance:")
print(rf_model$importance)

# list all classified tif files 
classified_files <- list.files(out_dir, pattern = "\\.tif$", full.names = TRUE)

# check if files exist
if(length(classified_files) == 0) {
  stop("No classified .tif files found in ", out_dir)
}

# -------------------------------------------------------------------
# calculating class percentages in each plot
# -------------------------------------------------------------------

# load shp first
plots <- st_read("~/Desktop/marshbirdsoutput/round_7/200m_buffer.shp") 

# create a lookup table with ID and Name before converting to terra
plot_lookup <- plots %>%
  st_drop_geometry() %>%  # remove geometry to create simple dataframe
  mutate(ID = row_number()) %>%  # create ID to match terra::extract output
  select(ID, Name)  # keep ID and Name columns

# convert to terra format for extract()
plots_terra <- vect(plots)

# initialize empty list to store results from each tile
all_extractions <- list()

# extract from each classified tile (memory efficient approach)
for(i in seq_along(classified_files)) {
  message("Extracting from tile ", i, " of ", length(classified_files), ": ", basename(classified_files[i]))
  
  tryCatch({
    # load one tile at a time
    tile_raster <- rast(classified_files[i])
    
    # extract values for plots that intersect this tile
    extracted_values <- terra::extract(tile_raster, plots_terra, df = TRUE)
    
    # rename the classification column to 'class' (it might be named differently)
    names(extracted_values)[2] <- "class"
    
    # only keep rows that actually have data (removes plots that don't intersect this tile)
    extracted_values <- extracted_values[!is.na(extracted_values$class), ]
    
    if(nrow(extracted_values) > 0) {
      all_extractions[[i]] <- extracted_values
    }
    
    # clean up memory
    rm(tile_raster)
    gc()
    
  }, error = function(e) {
    message("Error extracting from ", basename(classified_files[i]), ": ", e$message)
  })
}

# combine all extractions
all_values <- do.call(rbind, all_extractions)

# create class mapping
class_mapping <- data.frame(
  class_code = c(1, 2, 3, 4, 5, 6, 7, 8),  # numerical codes from raster
  class_name = c("hm", "lm", "ow", "ph", "rd", "up", "md", "sd")  # corresponding class names
)

# map numerical codes to class names

all_values <- all_values %>%
  mutate(class = as.numeric(class)) %>%  # convert factor to its underlying numeric level
  left_join(class_mapping, by = c("class" = "class_code")) %>%
  select(-class) %>%
  rename(class = class_name)

# calculate percent cover per plot
percents <- all_values %>%
  group_by(ID, class) %>%
  summarize(pixel_count = n(), .groups = "drop") %>%
  group_by(ID) %>%
  mutate(percent = pixel_count / sum(pixel_count)) %>%
  ungroup()

# join with plot names
percents_with_names <- percents %>%
  left_join(plot_lookup, by = "ID")

# pivot to wide format
percent_cover <- percents_with_names %>%
  select(ID, Name, class, percent) %>%
  pivot_wider(names_from = class, values_from = percent, values_fill = 0) %>%
  rename(plot_id = ID, plot_name = Name) %>%
  mutate(across(where(is.numeric), \(x) round(x, 2))) %>%
  relocate(plot_id, plot_name)  # put ID and Name columns first

# save to CSV
write.csv(percent_cover, "~/Desktop/marshbirdsoutput/round_7/percent_cover_by_plot.csv", row.names = FALSE)

