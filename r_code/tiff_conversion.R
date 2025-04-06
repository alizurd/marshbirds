Getting coordinates from raster

dem <- raster("usace2020_nj_ny_dem_J1158305_004_001.tif")
crs(dem)           # Coordinate Reference System
extent(dem)        # Bounding box (xmin, xmax, ymin, ymax)
res(dem)           # Pixel size
ncell(dem)         # Total number of pixels

coords <- coordinates(dem)        # matrix of x/y for each cell
values <- getValues(dem)          # elevation values
df <- data.frame(coords, elevation = values)

df <- as.data.frame(dem, xy = TRUE, na.rm = TRUE)
colnames(df) <- c("x", "y", "elevation")

View(df)
View(df[1:100, ])

Classifying the data

df <- as.data.frame(tiff, xy = TRUE, na.rm = TRUE)
colnames(df) <- c("x", "y", "elevation")

df$class <- with(df, ifelse([is.na](http://is.na/)(elevation), "nodata",
                            ifelse(elevation < 6.5, "Low Marsh",
                                   ifelse(elevation < 8, "High Marsh",
                                          ifelse(elevation < 10, "Road",
                                                 ifelse(elevation < 14, "Forest", "Other"))))))

final code

dem <- raster("usace2020_nj_ny_dem_J1158305_004_001.tif")
crs(dem)           # Coordinate Reference System
extent(dem)        # Bounding box (xmin, xmax, ymin, ymax)
res(dem)           # Pixel size
ncell(dem)         # Total number of pixels

coords <- coordinates(dem)        # matrix of x/y for each cell
values <- getValues(dem)          # elevation values
df <- data.frame(coords, elevation = values)

convert x and y, and elevation

# Assuming your current data frame is called `df`

# And columns are: x (easting), y (northing), elevation (in feet)

# Convert to sf object

df_sf <- st_as_sf(df, coords = c("x", "y"), crs = 6539)  # EPSG:6539 = NAD83 / New York Long Island (ftUS)

# Transform to WGS84 (latitude/longitude)

df_latlon <- st_transform(df_sf, crs = 4326)

# Extract lat/lon back to columns

df_coords <- df_latlon %>%
  mutate(lon = st_coordinates(.)[,1],
         lat = st_coordinates(.)[,2]) %>%
  st_drop_geometry()

df_coords <- df_coords %>%
  mutate(elevation_m = elevation * 0.3048)

df <- as.data.frame(dem, xy = TRUE, na.rm = TRUE)
colnames(df) <- c("elevation", "y", "x")

df <- df[, c("elevation", "y", "x")]  # if needed

View(df)
View(df_coords[1:100, ])