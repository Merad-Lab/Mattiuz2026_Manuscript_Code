library(tidyverse)
library(sf)

# Identify sample files
files_vec <- list.files("data/")[1:16]
`%notin%` <- function(x, y) !(`%in%`(x, y))

for (i in seq_along(files_vec)) {
  sampleID <- files_vec[i]
  print(paste("Processing cell density for:", sampleID))
  
  # Folder structure Setup
  output_dir <- paste0("celldensity/", sampleID, "/")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Load data
  sf_obj <- readRDS(paste0("preprocessing/", sampleID, "/", sampleID, "_clean.rds"))
  sf_obj <- subset(sf_obj, phenotype_clean %notin% c("Unknown", "Other"))
  
  sf_cells <- readRDS(paste0("preprocessing/", sampleID, "/", sampleID, "_cell_outline.rds"))
  sf_CCL19 <- readRDS(paste0("preprocessing/", sampleID, "/", sampleID, "_CCL19_geojson.rds"))
  sf_TLS <- readRDS(paste0("preprocessing/", sampleID, "/", sampleID, "_TLS_geojson.rds"))
  sf_tissue <- readRDS(paste0("preprocessing/", sampleID, "/", sampleID, "_tissue_outline.rds"))
  
  # Clean and merge CCL19 rows into cell data
  sf_CCL19_clean <- sf_CCL19 %>% 
    mutate(CellID = NA_real_, phenotype_clean = "CCL19_pos") %>% 
    select(CellID, phenotype_clean, geometry)
  
  # Classify CCL19 landmarks by TLS membership
  CCL19_with_TLS <- st_join(sf_CCL19_clean, sf_TLS[, c("id", "geometry")], join = st_intersects, left = TRUE)
  CCL19_with_TLS$TLS_status <- ifelse(is.na(CCL19_with_TLS$id), "non_TLS", "TLS")
  CCL19_with_TLS_df <- data.frame(phenotype_clean = CCL19_with_TLS$phenotype_clean, TLS_status = CCL19_with_TLS$TLS_status)
  
  # Compile consolidated metadata
  metadata_df <- sf_obj@meta.data %>% 
    dplyr::select(phenotype_clean, TLS_status) %>% 
    bind_rows(CCL19_with_TLS_df)
  
  # Density per area calculations
  tissue_area <- st_area(st_union(sf_tissue)) %>% as.numeric()
  TLS_area <- st_union(sf_TLS) %>% st_area() %>% as.numeric()
  nonTLS_area <- tissue_area - TLS_area
  
  # Calculate absolute frequencies and densities
  density_df <- metadata_df %>%
    group_by(TLS_status, phenotype_clean) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(
      area_pixels = ifelse(TLS_status == "TLS", TLS_area, nonTLS_area),
      area_mm2 = area_pixels / 1e6,
      cells_per_mm2 = count / area_mm2
    )
  
  # Export final clean summary dataset
  writexl::write_xlsx(density_df, path = paste0(output_dir, "differential_abundance_absolute_cell_density_mm2.xlsx"))
}
