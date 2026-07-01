library(tidyverse)
library(sf)

files_vec <- list.files("data/")[1:16]
`%notin%` <- function(x, y) !(`%in%`(x, y))

for (i in seq_along(files_vec)) {
  sampleID <- files_vec[i]
  print(paste("Processing neighborhood for:", sampleID))
  
  output_dir <- paste0("nhood/", sampleID, "/")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Load single-cell and regional matrices
  hd_obj <- readRDS(paste0("preprocessing/", sampleID, "/", sampleID, "_clean.rds"))
  cell_polygons <- readRDS(paste0("preprocessing/", sampleID, "/", sampleID, "_cell_outline.rds"))
  sf_TLS <- readRDS(paste0("preprocessing/", sampleID, "/", sampleID, "_TLS_geojson.rds"))
  sf_CCL19 <- readRDS(paste0("preprocessing/", sampleID, "/", sampleID, "_CCL19_geojson.rds"))
  
  # Set standardized identity mappings
  cell_polygons$barcode <- as.character(cell_polygons$CellID)
  
  sf_CCL19_clean <- sf_CCL19 %>%
    mutate(barcode = paste0("CCL19_", row_number()), phenotype_clean = "CCL19_pos") %>%
    select(barcode, phenotype_clean, geometry)
  
  # Establish TLS context mapping 
  CCL19_with_TLS <- st_join(sf_CCL19_clean, sf_TLS[, c("id", "geometry")], join = st_intersects, left = TRUE)
  CCL19_with_TLS_df <- data.frame(
    barcode = CCL19_with_TLS$barcode,
    TLS_status = ifelse(is.na(CCL19_with_TLS$id), "non_TLS", "TLS")
  )
  
  # Synthesize cellular centroids 
  cell_centroids <- bind_rows(cell_polygons %>% select(barcode, phenotype_clean, geometry), sf_CCL19_clean) %>% 
    filter(phenotype_clean %notin% c("Other", "Unknown")) %>% 
    st_centroid()
  
  # Assemble master metadata map
  meta_tls <- hd_obj@meta.data %>% 
    rownames_to_column("barcode") %>% 
    select(barcode, TLS_status) %>% 
    mutate(barcode = gsub(paste0(sampleID, "_"), "", barcode)) %>% 
    bind_rows(CCL19_with_TLS_df)
  
  # Query designated cell lineages for neighborhood evaluations
  ref_cell_types <- c("mature_cDC1", "mature_cDC2", "CCL19_pos")
  for (ref_type in ref_cell_types) {
    ref_cells <- cell_centroids %>% filter(phenotype_clean == ref_type)
    if (nrow(ref_cells) == 0) next
    
    # Calculate index matrices within radius bounds (30 units)
    pairs <- st_is_within_distance(ref_cells, cell_centroids, dist = 30)
    
    neighbor_table <- map2_dfr(seq_along(pairs), pairs, ~ {
      ref <- ref_cells[.x, ]
      neis <- cell_centroids[.y, ] %>% filter(barcode != ref$barcode)
      if (nrow(neis) == 0) return(NULL)
      
      tibble(
        ref_barcode = ref$barcode,
        ref_cluster = ref$phenotype_clean,
        neighbor_barcode = neis$barcode,
        neighbor_cluster = neis$phenotype_clean,
        distance_to_ref = as.numeric(st_distance(neis, ref))
      )
    })
    
    # Annotate spatial table with regional classification and save
    if (!is.null(neighbor_table) && nrow(neighbor_table) > 0) {
      neighbor_table <- neighbor_table %>% left_join(meta_tls, by = c("ref_barcode" = "barcode"))
      write.csv(neighbor_table, file = file.path(output_dir, paste0("neighbor_table_", ref_type, ".csv")), row.names = FALSE)
    }
  }
}
