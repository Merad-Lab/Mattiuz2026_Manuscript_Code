## Libraries ----
library(tidyverse)
library(semla)
library(scuttle)
library(Voyager)
library(SpatialFeatureExperiment)

############# RUN 1 ############# 
sampleID = "BIC18_1"
samples = "data/BIC18_1/filtered_feature_bc_matrix.h5"
imgs = "data/BIC18_1/spatial/tissue_hires_image.png"
spotfiles = "data/BIC18_1/spatial/tissue_positions_list.csv"
json = "data/BIC18_1/spatial/scalefactors_json.json"
infoTable = tibble(samples, imgs, spotfiles, json)

## Create object ----
se = ReadVisiumData(infoTable = infoTable)

## Add cell2location ----
cell2location_res = read.csv("data/BIC18_1/cell2location.txt")
colnames(cell2location_res) = gsub("q05cell_abundance_w_sf_", "", colnames(cell2location_res))
se@meta.data = cbind(se@meta.data, cell2location_res)

cell2loc_vec = se[[]]
cell2loc_vec = cell2loc_vec[, c(5:82)]

ref_celltypes = unique(colnames(cell2loc_vec))
ref_celltypes = ref_celltypes[1]
for (iter_celltype in ref_celltypes) {

  ## Voyager ----
  se@meta.data[, "array_row"] = se@tools$Staffli$x
  se@meta.data[, "array_col"] = se@tools$Staffli$y
  sfe = SpatialFeatureExperiment(list(counts = se@assays$Spatial$counts), colData = se@meta.data, spatialCoordsNames = c("array_row", "array_col"))
  sfe = sfe[rowSums(counts(sfe)) > 0, ]
  sfe = sfe[, colSums(counts(sfe)) > 0]
  sfe = logNormCounts(sfe)
  sfe = addPerFeatureQCMetrics(sfe)
  
  ## findSpatialNeighbors ----
  colGraph(sfe, "knn5_b") = findSpatialNeighbors(sfe, method = "knearneigh", dist_type = "idw", k = 5, style = "B")
  
  ## runUnivariate ---- 
  sfe@assays@data$logcounts["EGR2", ] = cell2loc_vec[, iter_celltype]
  sfe = runUnivariate(sfe, type = "localG_perm", features = "EGR2", colGraphName = "knn5_b", include_self = TRUE, swap_rownames = "symbol")
  lr = localResult(sfe, "localG_perm", "EGR2")
  se@meta.data[, paste0("Gi_", iter_celltype)] = lr$localG
  se_df = se[[]]
  
  gg1 = plotLocalResult(sfe, "localG_perm", features = "EGR2", colGeometryName = "centroids", divergent = TRUE, diverge_center = 0, size = 2)
  
  ## Visualize ----
  gg1 = ggplot(se_df, aes_string(x = "array_col", y = "array_row", color = paste0("Gi_", iter_celltype))) + geom_point(size = 2) + scale_color_gradientn(colours = pals::ocean.balance(100)) + ggdark::dark_theme_bw() + theme(legend.position = "bottom", title = element_blank()) + ggtitle("") 
  ggsave(gg1, filename = paste0("plots/", iter_celltype, "_hotspot_continous.svg"), width = 6, height = 8, create.dir = T)
}
