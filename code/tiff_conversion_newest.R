
# 1 - load and classify raster based on value ranges

install.packages("raster")
install.packages("tidyverse")
install.packages("sf")
install.packages("extractextractr")

library(raster)
library(sf)
library(dplyr)
library(exactextractr)

setwd("/Users/alyssabueno/Desktop/naip/tiffs")

# Load raster
r <- raster("bx_clip.tif")

# Classify into categories: low marsh = 1, high marsh = 2, open water = 3
r_classified <- calc(r, fun = function(x) {
  ifelse(x >= 0 & x <= 50, 1,
         ifelse(x > 50 & x <= 150, 2,
                ifelse(x > 150, 3, NA)))
})



# 2 - load lat/long points and buffer them 

# Load point data
points <- read.csv("latlong.csv")  # should have lat and lon columns
points_sf <- st_as_sf(points, coords = c("X", "Y"), crs = 32618)

View(points_sf)

# Reproject to match raster CRS
points_sf <- st_transform(points_sf, crs = crs(r))

# Create 200m buffer
buffers <- st_buffer(points_sf, dist = 200)



# 3 - calculate proportion of each class within each buffer

# Use exactextractr to get cell fractions in each buffer
results <- exact_extract(r_classified, buffers, fun = function(values, coverage_fractions) {
  total_area <- sum(coverage_fractions, na.rm = TRUE)
  low_marsh <- sum(coverage_fractions[values == 1], na.rm = TRUE) / total_area
  high_marsh <- sum(coverage_fractions[values == 2], na.rm = TRUE) / total_area
  open_water <- sum(coverage_fractions[values == 3], na.rm = TRUE) / total_area
  return(data.frame(low_marsh = low_marsh, high_marsh = high_marsh, open_water = open_water))
})



# 4 - combine with coordinates and export 

# Get lat/lon back
coords <- st_coordinates(st_centroid(buffers))

# Combine with results
final_df <- cbind(data.frame(x = coords[,2], y = coords[,1]), results)

# Save as CSV
write.csv(final_df, "classified_summary_bx.csv", row.names = FALSE)

