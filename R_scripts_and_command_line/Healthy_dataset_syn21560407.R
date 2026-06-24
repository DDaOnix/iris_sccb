# Healty dataset from https://www.synapse.org/Synapse:syn21560407
library(Seurat)
library(reticulate)
library(anndata)
library("Nebulosa")
library("SCpubr")

# Convert the .h5as into seurat object
data <- read_h5ad("/Volumes/prj_id_iris/IRIS/Healthy_ref_GSE143868/droplet_normal_lung_blood_scanpy.20200205.RC4.h5ad")
data <- CreateSeuratObject(counts = t(as.matrix(data$X)), meta.data = data$obs,min.features = 500, min.cells = 30)
saveRDS(data,"/Volumes/prj_id_iris/IRIS/Healthy_ref_GSE143868/droplet_normal_lung_blood_scanpy.20200205.RC4.rds")

# QC and filtering
data$mito.percent <- PercentageFeatureSet(data, pattern = '^MT-')
View(data@meta.data)
# explore QC

# filter
data.filtered <- subset(data, subset = nCount_RNA > 800 &
                                   nFeature_RNA > 200 & 
                                   nFeature_RNA < 5000 & # Set an upper limit of genes expressed per cell as it can limit the inclusion of doublets or more-blets
                                   mito.percent < 5) # <10 suggested by Osorio D, Cai JJ. 

# Systematic determination of the mitochondrial proportion in human 
# and mice tissues for single-cell RNA-sequencing data quality control. 
# Bioinformatics. 2021 May 17;37(7):963-967. 
# doi: 10.1093/bioinformatics/btaa751. PMID: 32840568; PMCID: PMC8599307.

# standard workflow steps
data.filtered <- NormalizeData(data.filtered)
data.filtered <- FindVariableFeatures(data.filtered)
data.filtered <- ScaleData(data.filtered)
data.filtered <- RunPCA(data.filtered)
ElbowPlot(data.filtered)
data.filtered <- RunUMAP(data.filtered, dims = 1:20, reduction = 'pca')

data.filtered <- saveRDS(data.filtered, "/Volumes/prj_id_iris/IRIS/Healthy_ref_GSE143868/data.filtered_healthyall.rds")
data.filtered <- readRDS("/Volumes/prj_id_iris/IRIS/Healthy_ref_GSE143868/data.filtered_healthyall.rds")

plot_density(data.filtered, c("IL5RA", 'CSF2RB'), reduction = "umap")

plot_density(data.filtered, c("IGF2BP3"), reduction = "umap")

DimPlot(data.filtered, reduction = 'umap', group.by = 'free_annotation', label = TRUE, raster=FALSE) + NoLegend()
DimPlot(data.filtered, reduction = 'umap', group.by = 'magnetic.selection', label = TRUE, raster=FALSE) 

RidgePlot(data.filtered, features = c("IGF2BP3"), group.by = 'magnetic.selection', ncol = 2)
VlnPlot(data.filtered, features = c("IGF2BP3"), group.by = 'magnetic.selection', pt.size = 0)
VlnPlot(data.filtered, features = c("IGF2BP3"), group.by = 'free_annotation', pt.size = 0) + NoLegend
FeaturePlot(data.filtered, features = c("IGF2BP3"), min.cutoff = 'q9', order = TRUE)

DotPlot(data.filtered, features = c("IL5RA", 'CSF2RB'), group.by = 'free_annotation') + RotatedAxis() +
  theme(axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=7),
        panel.grid.major.y = element_line(color = "grey", size = 0.2))


# Create a new column in meta.data without the "_P1", "P2" or "_P3" suffix
data.filtered@meta.data$cell.type <- gsub("_P[123]$", "", data.filtered@meta.data$free_annotation)
View(data.filtered@meta.data)
DimPlot(data.filtered, reduction = 'umap', group.by = 'cell.type', label = TRUE, raster=FALSE, repel = TRUE) + NoLegend()
