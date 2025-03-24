library(Seurat)
library(spacexr)
library(Matrix)
library(dplyr)

## Loading the single cell Reference

path_raw_sc <- 'your_inpath_folder' # Replace with the location of your single-nucleus RNA sequencing data

expression_matrix_sc <- ReadMtx(
  mtx = paste0(path_raw_sc, 'matrix.mtx'), features = paste0(path_raw_sc, 'genes.tsv'),
  cells = paste0(path_raw_sc, 'barcodes.tsv')
)

# Create a Seurat object
seurat_raw_sc <- CreateSeuratObject(counts = expression_matrix_sc)

# I need the cell type annotations which are on a anndata object. 
# I will export in python and import here the barcodes together with the cell
# types

annotations <- 
  read.csv('CellTypeAnnotations_RCTD.csv')

## We remove mixed and BT cell doublets

exclude_celltypes <- c('B T cell doublet', 'mixed')
annotations_final <- annotations %>% dplyr::filter(!(celltype_merged %in% exclude_celltypes))

seurat_subset_sc <- subset(seurat_raw_sc, cells = annotations_final$CELL)

seurat_subset_sc@meta.data$CELL <- colnames(seurat_subset_sc)

seurat_subset_sc@meta.data <-  
  seurat_subset_sc@meta.data %>% dplyr::inner_join(annotations, by ='CELL')

cell_types <- as.factor(seurat_subset_sc@meta.data$celltype_merged)
names(cell_types) <-  seurat_subset_sc$CELL

nUMI_sc <- seurat_subset_sc@meta.data$nCount_RNA
names(nUMI_sc) <- seurat_subset_sc$CELL

### Create the Reference object
reference <- Reference(counts = GetAssayData(object = seurat_subset_sc, assay = "RNA", layer = "counts"), 
                       cell_types = cell_types, nUMI = nUMI_sc)


rm(annotations, annotations_final, expression_matrix_sc, seurat_raw_sc, seurat_subset_sc)

# After removing objects, you can call garbage collection explicitly
gc()

### Reading the spatial data. I will need to filter spots that did not pass the quality control

path_raw_st <- 'your_inpath_folder' # Replace with the location of your spatial transcriptomics data'

samples <- list.files(path = path_raw_st)

seurat_sp_list <- list()
resultsdir <- '/deconvolution/RCTD/'
dir.create(resultsdir)

print(paste0("Number of cores :", parallel::detectCores()))


for (current_sample in samples){ 
  print(current_sample)
  
  current_seurat_object <- Load10X_Spatial(paste0(path_raw_st, current_sample))
  
  coords_sp <- GetTissueCoordinates(current_seurat_object)
  
  nUMI_sp <- current_seurat_object@meta.data$nCount_Spatial
  names(nUMI_sp) <- colnames(current_seurat_object)
  
  puck <-SpatialRNA(coords_sp, 
                    counts=GetAssayData(object = current_seurat_object, 
                    assay = "Spatial", 
                    layer = "counts"), 
                    nUMI_sp)
  
  myRCTD <- create.RCTD(puck, reference, 
                        max_cores = parallel::detectCores(), 
                        CELL_MIN_INSTANCE = 10)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
  
  results <- myRCTD@results
  norm_weights = normalize_weights(results$weights) 
  cell_type_names <- myRCTD@cell_type_info$info[[2]]
  spatialRNA <- myRCTD@spatialRNA
  
  ## Adding the norm weights to the seurat object. 
  current_seurat_object_subset <- 
    subset(current_seurat_object, cells = rownames(norm_weights))

  current_seurat_object_subset <- 
    AddMetaData(object = current_seurat_object_subset, metadata = norm_weights)
  
  seurat_sp_list[[current_sample]] <- current_seurat_object_subset
  
}

saveRDS(seurat_sp_list, paste0(resultsdir,"Seurat_object_RCTD_results.rds"))


resultsdir <- '/deconvolution/RCTD/'
seurat_sp_list <- readRDS(paste0(resultsdir,"Seurat_object_RCTD_results.rds"))


  
  


