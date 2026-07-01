library(tidyverse)
library(sf)

files_vec <- list.files("data/")[1:16]
`%notin%` <- function(x, y) !(`%in%`(x, y))

for (i in seq_along(files_vec)) {
  sampleID <- files_vec[i]
  print(paste("Processing distances for:", sampleID))
  
  output_dir <- paste0("distance/", sampleID, "/")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Load spatial layers
  sf_cells <- readRDS(paste0("preprocessing/", sampleID, "/", sampleID, "_cell_outline.rds"))
  sf_cells <- subset(sf_cells, phenotype_clean %notin% c("Unknown", "Other"))
  sf_CCL19 <- readRDS(paste0("preprocessing/", sampleID, "/", sampleID, "_CCL19_geojson.rds"))
  
  # Compute absolute distance to the nearest target landmark feature
  nearest_idx <- st_nearest_feature(sf_cells, sf_CCL19)
  distances <- st_distance(sf_cells, sf_CCL19[nearest_idx, ], by_element = TRUE)
  
  sf_cells <- sf_cells %>% 
    mutate(
      dist_to_CCL19 = as.numeric(distances),
      proximity = case_when(dist_to_CCL19 <= 50 ~ "proximal", dist_to_CCL19 > 50 ~ "distal")
    )
  
  # Segment distances into even breaks (bins)
  breaks_vec <- seq(0, max(sf_cells$dist_to_CCL19, na.rm = TRUE), length.out = 11)
  sf_cells <- sf_cells %>% 
    mutate(dist_bin = cut(dist_to_CCL19, breaks = breaks_vec, include.lowest = TRUE, dig.lab = 5))
  
  # Drop geometry prior to heavy summary calculations
  cell_summary_bins <- sf_cells %>%
    st_drop_geometry() %>% 
    group_by(phenotype_clean, dist_bin) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(dist_bin) %>%
    mutate(percent = count / sum(count) * 100) %>%
    ungroup()
  
  # Export core structural datasets
  saveRDS(sf_cells, file = paste0(output_dir, sampleID, "_distance_to_CCL19_sf_cells.rds"))
  saveRDS(cell_summary_bins, file = paste0(output_dir, sampleID, "_distance_summary_bins.rds"))
}
