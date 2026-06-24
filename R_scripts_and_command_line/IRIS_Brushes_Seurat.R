# IRIS BB pre-post
# script to integrate scRNA-Seq datasets to correct for batch effects
setwd("~/Desktop/scRNA-seq_Brush1/IRIS/Filtered_feature_bc_matrix/")


# load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)

# get filtered_feature_bc_matrix location
dirs <- list.dirs(path = './', recursive = F, full.names = F)

for(x in dirs){
  name <- gsub('_filtered_feature_bc_matrix','', x)
  
  cts <- ReadMtx(mtx = paste0('./',x,'/matrix.mtx.gz'),
                 features = paste0('./',x,'/features.tsv.gz'),
                 cells = paste0('./',x,'/barcodes.tsv.gz'))
  
  # create seurat objects
  assign(name, CreateSeuratObject(counts = cts))
}




# merge datasets
ls()
merged_seurat <- merge(P11_A, y = c(P11_B, P13_A, P13_B, P16_A, P16_B, P18_A, P18_B, P19_A, P19_B, P2_A, P2_B, P20_A, P20_B,
                                    P3_A, P3_B, P7_A, P7_B, P8_A, P8_B, P9_A, P9_B),
                       add.cell.ids = ls()[4:25], #this is selecting the elements output from ls() that is each seurat object created above
                       project = 'IRIS_Br')


merged_seurat

# QC & filtering -----------------------

View(merged_seurat@meta.data)
# create a sample column
merged_seurat$sample <- rownames(merged_seurat@meta.data)

# split sample column
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('Patient', 'Type', 'Barcode'), 
                                    sep = '_')

# Check all patient have been merged 
unique(merged_seurat@meta.data$Patient)
# Check all conditions have been mergrd
unique(merged_seurat@meta.data$Type)

# calculate mitochondrial percentage
merged_seurat$mitoPercent <- PercentageFeatureSet(merged_seurat, pattern='^MT-')

# explore QC
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "^MT-")
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Save the object
saveRDS(merged_seurat, "~/Desktop/scRNA-seq_Brush1/IRIS/merged_seurat_obj.rds")



# filtering
merged_seurat_filtered <- subset(merged_seurat, subset = nCount_RNA > 800 &
                                   nFeature_RNA > 500 &
                                   mitoPercent < 10)

merged_seurat_filtered

merged_seurat




# perform standard workflow steps to figure out if we see any batch effects --------
merged_seurat_filtered <- NormalizeData(object = merged_seurat_filtered)
merged_seurat_filtered <- FindVariableFeatures(object = merged_seurat_filtered)
merged_seurat_filtered <- ScaleData(object = merged_seurat_filtered)
merged_seurat_filtered <- RunPCA(object = merged_seurat_filtered)
ElbowPlot(merged_seurat_filtered)
merged_seurat_filtered <- FindNeighbors(object = merged_seurat_filtered, dims = 1:20)
merged_seurat_filtered <- FindClusters(object = merged_seurat_filtered)
merged_seurat_filtered <- RunUMAP(object = merged_seurat_filtered, dims = 1:20)

# Save the object
saveRDS(merged_seurat_filtered, "~/Desktop/scRNA-seq_Brush1/IRIS/merged_seurat_filtered_obj.rds")
# plot
p1 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'Patient')
p2 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'Type',
              cols = c('grey','royalblue'))

grid.arrange(p1, p2, ncol = 2, nrow = 2)

DimPlot(merged_seurat_filtered, reduction = 'umap')

# perform integration to correct for batch effects ------
# obj.list <- SplitObject(merged_seurat_filtered, split.by = 'Patient')
# for(i in 1:length(obj.list)){
#   obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
#   obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
# }
############
obj.list <- SplitObject(merged_seurat_filtered, split.by = 'Patient')
for(i in 1:length(obj.list)){
  obj.list[[i]] <- JoinLayers(object = obj.list[[i]]) #this line changes from the code above. This solves and error (GetAssayData doesn't work for multiple layers in v5 assay)
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}
#############



# select integration features
features <- SelectIntegrationFeatures(object.list = obj.list)

# find integration anchors (CCA)
anchors <- FindIntegrationAnchors(object.list = obj.list,
                                  anchor.features = features)

# integrate data
seurat.integrated <- IntegrateData(anchorset = anchors)


# Scale data, run PCA and UMAP and visualize integrated data
seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated)
seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:50)


p3 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'Patient')
p4 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'Type',
              cols = c('grey','royalblue'))


grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)

















