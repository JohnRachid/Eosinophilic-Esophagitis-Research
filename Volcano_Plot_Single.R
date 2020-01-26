setwd('C:\\Users\\John\\Documents\\EOE Project\\')
getwd()

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("multtest")
BiocManager::install("DESeq2")
BiocManager::install("EnhancedVolcano")
BiocManager::install("airway")
install.packages('Seurat')
devtools::install_github("kevinblighe/EnhancedVolcano")


library(Seurat)
library(EnhancedVolcano)
library(airway)
library(magrittr)
library(DESeq2)


lung <- gzfile("lung.rds")
lung_RDS<- readRDS(lung)
#begin csv of mast cells 
lung_meta <- lung_RDS@meta.data
lung_meta_df <- as.data.frame(lung_meta)
lung_df_mast <- lung_meta_df[lung_meta_df$Celltypes == "Mast_cells",] #here we have all of the mast cell values. However we need to get a col for p value and logfolds
write.csv(lung_df_mast, "lung_mast_csv.csv")
#end csv of mast cells

colnames(lung_RDS)
lung_RDS[["RNA"]]@counts



# pbmc <- CreateSeuratObject(pbmc,  project = "SeuratProject",assay = "RNA",)
#pbmc <- NormalizeData(lung_RDS, normalization.method = "LogNormalize", scale.factor = 10000)
#lung.markers <- FindAllMarkers(lung_RDS,assay = "RNA",) #goes based on columns
GetAssayData(object = lung_RDS, slot = 'scale.data')[1:3, 1:3]
lung.markers <- FindAllMarkers(lung_RDS,assay = "RNA") #todo runs out of ram


lung_plot <- EnhancedVolcano(lung.markers,
                                   lab = rownames(lung.markers),
                                   x = 'avg_logFC',
                                   y = 'p_val',
                                   xlim = c(-3.5, 3.5),
                                   title= "FindAllMarkers Lung")

lung_plot


lung_plot_boxed <- EnhancedVolcano(lung.markers,
                                         lab = rownames(lung.markers),
                                         x = 'avg_logFC',
                                         y = 'p_val',
                                         xlim = c(-3.5, 3.5),
                                         boxedLabels = TRUE,
                                         title= "FindAllMarkers Lung Boxed")

lung_plot_boxed




oesophagus <- gzfile("oesophagus.rds")
oesophagus_RDS <- readRDS(oesophagus)
#begin create csv of mast cells
oesophagus_meta <- oesophagus_RDS@meta.data
oesophagus_meta_df <- as.data.frame(oesophagus_meta)
oesophagus_df_mast <- oesophagus_meta_df[oesophagus_meta_df$Celltypes == "Mast_cell",]
write.csv(oesophagus_df_mast, "oesophagus_mast_csv.csv")
#end create csv of mast cells
oesophagus.markers <- FindAllMarkers(oesophagus_RDS,assay = "RNA",) #todo runs out of ram





oesophagus_plot <- EnhancedVolcano(oesophagus.markers,
                lab = rownames(oesophagus.markers),
                x = 'avg_logFC',
                y = 'p_val',
                xlim = c(-3.5, 3.5),
                title= "FindAllMarkers oesophagus")

oesophagus_plot


oesophagus_plot_boxed <- EnhancedVolcano(oesophagus.markers,
                                   lab = rownames(oesophagus.markers),
                                   x = 'avg_logFC',
                                   y = 'p_val',
                                   xlim = c(-3.5, 3.5),
                                   boxedLabels = TRUE,
                                   title= "FindAllMarkers Oesophagus Boxed")

oesophagus_plot_boxed