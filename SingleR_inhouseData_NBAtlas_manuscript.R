# Inputs ------------------------------------------------------------------
test_dataset_path <- "/scratch/gent/vo/000/gvo00027/projects/Single_Cell_Neuroblastoma/NBAtlas/05f_SingleR_inhouseData_NBAtlas/Robjects/05f_SingleR_inhouseData_input_seurat_NBAtlas.rds"# path to seurat test obj (needs to be log-normalized)

ref_dataset_path <- "/scratch/gent/vo/000/gvo00027/projects/Single_Cell_Neuroblastoma/NBAtlas/03c_post_scVI_R_plots/Robjects/seuratObj_NBAtlas_share_20230525.rds"# path to seurat ref obj (with annotation, needs to be log-normalized)
ref_annot_column <- "annot_NBN_scarches"# reference annotation column to use for classification

output_folder <- "05f_SingleR_inhouseData_NBAtlas"
project <- "NBAtlas"

output_name <- "inhouseData_vs_NBAtlas_" #part of Robject output name

nCores <- 18

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

ref <- readRDS(ref_dataset_path)
ref$reference_annotation <- ref@meta.data[ref_annot_column][,1]

message("Starting SingleR...")
SingleR_pred <- SingleR(test = test@assays$RNA@data,  assay.type.test = "logcounts", ref = ref@assays$RNA@data, assay.type.ref = "logcounts", labels = ref$reference_annotation, BPPARAM=MulticoreParam(nCores))

saveRDS(object = SingleR_pred, file = paste0(output_Robjects, output_name, "SingleR_predictions", project, ".rds"))