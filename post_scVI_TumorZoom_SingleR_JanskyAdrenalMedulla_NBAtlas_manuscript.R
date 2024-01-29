
# Inputs ------------------------------------------------------------------
test_dataset_path <- "06b_post_scVI_TumorZoom_R_import_NBAtlas/Robjects/06b_post_scVI_TumorZoom_R_import_TumorZoom_scVI_covSample_seuratObj_NormOnly_NBAtlas.rds" # path to seurat test obj (needs to be log-normalized)

ref_dataset_path <- "/scratch/gent/vo/000/gvo00027/data/Single_Cell_Neuroblastoma/NB_human/Author_processed/Jansky2021_adrenal_medulla/adrenal_medulla_Seurat.RDS" # path to seurat ref obj (with annotation, needs to be log-normalized)
#ref_annot_column <- ""# reference annotation column to use for classification

output_folder <- "06c_post_scVI_TumorZoom_R_plots_NBAtlas/"
step <- "06c_post_scVI_TumorZoom_R_plots_covSample"
project <- "NBAtlas"

output_name <- "Jansky2021_AdrenalMedulla" #part of Robject output name

nCores <- 8

# Libs --------------------------------------------------------------------
suppressPackageStartupMessages({
library(Seurat)
library(SingleR)
library(BiocParallel)
})
  
set.seed(123)

# Run ---------------------------------------------------------------------
output_Robjects <- paste0(output_folder, "/Robjects/")
dir.create(output_Robjects)

test <- readRDS(test_dataset_path)
#test_sce <- as.SingleCellExperiment(test) #GetAssayData(test,slot = "data")

#rm(test)
#gc()

ref <- readRDS(ref_dataset_path)

ref$reference_annotation <- Idents(ref) # only in this case

SingleR_pred <- SingleR(test = test@assays$RNA@data,  assay.type.test = "logcounts", ref = ref@assays$RNA@data, assay.type.ref = "logcounts", labels = ref$reference_annotation, BPPARAM=MulticoreParam(nCores))

saveRDS(object = SingleR_pred, file = paste0(output_Robjects, step, "_", output_name, "SingleR_predictions", project, ".rds"))