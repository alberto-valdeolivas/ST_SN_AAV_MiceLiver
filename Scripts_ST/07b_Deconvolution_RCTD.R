library(Seurat)
# library(Matrix)
# library(dplyr)

resultsdir <- '/deconvolution/RCTD/'
seurat_sp_list <- readRDS(paste0(resultsdir,"Seurat_object_RCTD_results.rds"))


seurat_sp_list[[1]]@meta.data


for (i in seq_along(seurat_sp_list)) {
  # Extract the current data frame
  current_df <- seurat_sp_list[[i]]@meta.data
  
  current_df$orig.ident <- names(seurat_sp_list)[i]
  current_df$spot_id <- rownames(current_df)
  current_df$spot_id_barcode <- paste(current_df$spot_id, current_df$orig.ident, sep = "-")
  
  # Define the file name dynamically
  file_name <- paste0(resultsdir, "deconv_RCTD_results_", names(seurat_sp_list)[i], ".csv")
  
  # Write the extracted data to a .csv file
  write.csv(current_df, file = file_name, row.names = FALSE)
} 
