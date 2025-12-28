#!/usr/bin/env Rscript
# Calculate percent cover by vegetation class per plot and write a final CSV

suppressPackageStartupMessages({
  library(optparse)
  library(terra)
  library(sf)
  library(dplyr)
  library(tidyr)
  library(readr)
})

option_list <- list(
  make_option(c("-d", "--classified_dir"), type = "character", default = "~/Desktop/marshbirdsoutput",
              help = "Directory containing classified .tif files (default: %default)", metavar = "path"),
  make_option(c("-f", "--files"), type = "character", default = NULL,
              help = "Comma-separated list of classified .tif files to process (overrides --classified_dir)", metavar = "files"),
  make_option(c("-s", "--plots_shp"), type = "character", default = "~/Desktop/r_inputs/200m_buffer_clip.shp",
              help = "Plot shapefile with buffers (default: %default)", metavar = "path"),
  make_option(c("-o", "--output"), type = "character", default = "data/percent_cover_by_plot.csv",
              help = "Output CSV path (default: %default)", metavar = "path")
)

opts <- parse_args(OptionParser(option_list = option_list))

classified_files <- if (!is.null(opts$files)) {
  strsplit(opts$files, ",")[[1]] |> trimws()
} else {
  d <- path.expand(opts$classified_dir)
  if (!dir.exists(d)) d <- opts$classified_dir
  list.files(d, pattern = "_classified\\.tif$|classified_output\\.tif$|_model.*_classified\\.tif$|\\.tif$",
             full.names = TRUE)
}

if (length(classified_files) == 0) {
  stop("No classified .tif files found. Provide --files or point --classified_dir to files.")
}

plots_shp <- path.expand(opts$plots_shp)
if (!file.exists(plots_shp)) stop("Plots shapefile not found: ", plots_shp)

plots <- sf::st_read(plots_shp, quiet = TRUE)
plots_vect <- terra::vect(plots)

# default mapping for class codes -> class names (adjust if your codes differ)
class_labels <- tibble::tibble(
  class = c(1, 2, 3, 4, 5, 6, 7),
  class_name = c("hm", "lm", "ow", "ph", "rd", "up", "md")
)

all_percent_cover <- list()

for (f in classified_files) {
  message("Processing: ", f)
  r <- tryCatch(terra::rast(f), error = function(e) NULL)
  if (is.null(r)) {
    warning("Could not read raster: ", f); next
  }

  vals <- terra::extract(r, plots_vect, df = TRUE)
  if (nrow(vals) == 0) {
    warning("No extracted values for: ", f); next
  }

  # raster value column typically is the last column (ID is first)
  val_col <- setdiff(names(vals), "ID")[1]
  if (is.null(val_col)) {
    warning("Could not find raster value column in extracted table for: ", f); next
  }

  vals <- vals %>%
    rename(class = !!sym(val_col)) %>%
    filter(!is.na(class))

  percents <- vals %>%
    group_by(ID, class) %>%
    summarize(pixel_count = n(), .groups = "drop") %>%
    group_by(ID) %>%
    mutate(percent = pixel_count / sum(pixel_count)) %>%
    ungroup()

  # join names when available
  percents <- percents %>%
    mutate(class = as.numeric(class)) %>%
    left_join(class_labels, by = "class") %>%
    mutate(class = ifelse(is.na(class_name), as.character(class), class_name)) %>%
    select(ID, class, percent)

  percent_cover <- percents %>%
    pivot_wider(names_from = class, values_from = percent, values_fill = 0) %>%
    rename(plot_id = ID) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
    mutate(source = basename(f))

  all_percent_cover[[f]] <- percent_cover
}

if (length(all_percent_cover) == 0) stop("No percent cover tables generated.")

final <- bind_rows(all_percent_cover) %>%
  arrange(source, plot_id)

# ensure output directory exists
out_path <- opts$output
out_dir <- dirname(out_path)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

readr::write_csv(final, out_path)
message("Wrote percent cover CSV to: ", out_path)

# end
