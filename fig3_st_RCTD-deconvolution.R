####################################
# RCTD for 10x visium spot deconvolution
# Dataset: D15
# Yiming 2022/6/21
####################################

library(Seurat)
library(spacexr)
library(Matrix)
library(tidyverse)
library(dplyr)


setwd("/usersdata/yiming/VIO/Ranalysis/RCTD_june")

# ----- load spatial and sc data -----
sample_ST <- Load10X_Spatial(data.dir = "/usersdata/yiming/VIO/spaceranger_run/D1/outs",
                             filename = "filtered_feature_bc_matrix.h5", assay = "Spatial",
                             slice = "slice4", filter.matrix = TRUE, to.upper = FALSE, image = NULL)
sample_sc = readRDS("/usersdata/yiming/VIO/Ranalysis/apr2/singlet_mnn_anno.rds")


# ----- prepare sc reference -----
counts <- GetAssayData(sample_sc, slot = "counts")
meta_data <- sample_sc@meta.data 
cell_types <- meta_data$celltype; names(cell_types) <- rownames(meta_data) 
cell_types <- as.factor(cell_types) 
nUMI <- meta_data$nCount_RNA; names(nUMI) <- rownames(meta_data) 
reference_sc <- Reference(counts, cell_types, nUMI)


# ----- check ref & save -----
print(dim(reference_sc@counts))
table(reference_sc@cell_types)
saveRDS(reference_sc, 'reference_sc.rds')


# ----- create spatial object -----
coords <- sample_ST@images[["slice4"]]@coordinates[, c('row', 'col')]
counts <- GetAssayData(sample_ST, slot = "counts")
nUMI <- colSums(counts)
sample_ST <- SpatialRNA(coords, counts, nUMI)


# ----- run RCTD -----
reference_sc <- readRDS('/usersdata/yiming/VIO/Ranalysis/RCTD_june/reference_sc.rds')
# default parameters
myRCTD <- create.RCTD(sample_ST, reference_sc, max_cores = 4, 
                      gene_cutoff = 0.000125, fc_cutoff = 0.5,
                      gene_cutoff_reg = 2e-04, fc_cutoff_reg = 0.75)
# run 
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')


# ----- check RCTD results -----

results <- myRCTD@results
norm_weights = normalize_weights(results$weights) 
cell_type_names <- myRCTD@cell_type_info$info[[2]] 
sample_ST <- myRCTD@spatialRNA


# ----- check cell type proportion-----

# select the top 2 weighted value for each spot
norm_weights
df_norm_wei <- as.data.frame(norm_weights)
write.csv(df_norm_wei, 'slice4_norm_wei.csv') 

