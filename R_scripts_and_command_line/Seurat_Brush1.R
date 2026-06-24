# SEURAT scRNA-seq data analysis guide from https://satijalab.org/seurat/articles/brush13k_tutorial
library(dplyr)
library(Seurat)
library(patchwork)

# Load the Brush1 dataset
brush1.data <- Read10X(data.dir = "~/Desktop/scRNA-seq_Brush1/filtered_feature_bc_matrix/")
# Initialize the Seurat object with the raw (non-normalized data).
brush1 <- CreateSeuratObject(counts = brush1.data, project = "brush1", min.cells = 3, min.features = 200)
brush1

# Standard pre-processing workflow
# QC and selecting cells for further analysis
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
brush1[["percent.mt"]] <- PercentageFeatureSet(brush1, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(brush1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(brush1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(brush1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
brush1 <- subset(brush1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalizing the data
brush1 <- NormalizeData(brush1, normalization.method = "LogNormalize", scale.factor = 10000)
# For clarity, in this previous line of code (and in future commands), we provide the default values for certain parameters in the function call. 
# However, this isn’t required and the same behavior can be achieved with:
#  brush1 <- NormalizeData(brush1)

# Identification of highly variable features (feature selection)
# We next calculate a subset of features that exhibit high cell-to-cell variation in the dataset (i.e, they are highly expressed in some cells, and lowly expressed in others). 
brush1 <- FindVariableFeatures(brush1, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(brush1), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(brush1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#  plot1 + plot2 - gives me error because of window size
plot1
plot2

# Scaling the data
# Next, we apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA. 
all.genes <- rownames(brush1)
brush1 <- ScaleData(brush1, features = all.genes)

# Perform linear dimensional reduction
# Next we perform PCA on the scaled data. 
# For the first principal components, Seurat outputs a list of genes with the most positive and negative loadings, 
# representing modules of genes that exhibit either correlation (or anti-correlation) across single-cells in the dataset.
brush1 <- RunPCA(brush1, features = VariableFeatures(object = brush1))

VizDimLoadings(brush1, dims = 1:2, reduction = "pca")

DimPlot(brush1, reduction = "pca") + NoLegend()

DimHeatmap(brush1, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(brush1, dims = 1:15, cells = 500, balanced = TRUE)

# Determine the ‘dimensionality’ of the dataset
# To overcome the extensive technical noise in any single feature for scRNA-seq data, 
# Seurat clusters cells based on their PCA scores, with each PC essentially representing 
# a ‘metafeature’ that combines information across a correlated feature set. 
# The top principal components therefore represent a robust compression of the dataset. 
# However, how many components should we choose to include? 10? 20? 100?
ElbowPlot(brush1)

# Cluster the cells
brush1 <- FindNeighbors(brush1, dims = 1:10)
brush1 <- FindClusters(brush1, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(brush1), 5)

# Run non-linear dimensional reduction (UMAP/tSNE)
brush1 <- RunUMAP(brush1, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(brush1, reduction = "umap", label = TRUE)
# Save the object at this point so that it can easily be loaded back in without having to rerun the computationally intensive steps performed above, 
# or easily shared with collaborators.
############## saveRDS(brush1, file = "~/Desktop/scRNA-seq_Brush1/brush1_obj.rds")

# To load it
brush1 <- LoadSeuratRds("~/Desktop/scRNA-seq_Brush1/brush1_obj.rds")


# Finding differentially expressed features (cluster biomarkers)
# Seurat can help you find markers that define clusters via differential expression (DE).
# FindAllMarkers() identifies positive and negative markers for all clusters, but you can also test groups of clusters vs. each other, or against all cells.
# find markers for every cluster compared to all remaining cells, report only the positive ones
clusters.markers <- FindAllMarkers(brush1, only.pos = TRUE)
clusters.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
  
###### THIS IS FOR SINGLE CLUSTER ONLY - IF YOU RUN IT< IT WILL DISRUPT THE GENERAL HEATMAP - JUMP TO "CONTINUE HERE"
# For single cluster biomarkers:
# find all markers of cluster 2
#   cluster0.markers <- FindMarkers(brush1, ident.1 = 0)
#   head(cluster0.markers, n = 20)
# Seurat has several tests for differential expression which can be set with the test.use parameter (see our DE vignette for details). 
# For example, the ROC test returns the ‘classification power’ for any individual marker (ranging from 0 - random, to 1 - perfect).
clusters.markers <- FindMarkers(brush1, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(brush1, features = c("TPSAB1", "TPSB2", "CPA3", "GATA2", "KIT", "MS4A2", "ENO2", "HDC", "HPGDS", "RGS13")) # mast cells markers

FeaturePlot(brush1, features = c("TPSAB1", "TPSB2", "CPA3", "GATA2", "KIT", "MS4A2", "ENO2", "HDC", "HPGDS", "RGS13")) # mast cells markers

# What I call Macrophages
FeaturePlot(brush1, features = c("MEF2C", "HLA-DQA1", "SLC8A1", "AIF2", "HLA-DQB1", "LYX", "MRC1", "MS4A7", "IFI30", "C1QB"))
##
########

## CONTINUE HERE

clusters.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(brush1, features = top10$gene) + NoLegend()

# Assigning cell type identity to clusters
new.cluster.ids <- c("Immune cells", "Immune cells", "Mast cells", "Ciliated", "NA", "Macrophages", "Mast cells", "Basal")
names(new.cluster.ids) <- levels(brush1)
brush1 <- RenameIdents(brush1, new.cluster.ids)
DimPlot(brush1, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# To calculate the number of cells per cluster from https://satijalab.org/seurat/archive/v3.0/interaction_vignette
table(Idents(brush1))
# What proportion of cells are in each cluster?
prop.table(table(Idents(brush1)))
# Plot UMAP, coloring cells by cell type (currently stored in object@ident)
DimPlot(brush1, reduction = "umap")

# Dot plots - the size of the dot corresponds to the percentage of cells expressing the
# feature in each cluster. The color represents the average expression level
DotPlot(brush1, features = c("TPSAB1", "TPSB2", "CPA3", "GATA2", "KIT", "MS4A2", "ENO2", "HDC", "HPGDS", "RGS13")) + RotatedAxis() # mast cells markers




### 
Idents(brush1)
# 
FeaturePlot(brush1, features = c("POUF2F3", "TRPM5", "PSMD5", "NGF"))







## Try to detect doublets with DoubletFinder 
# remotes::install_github('chris-mcginnis-ucsf/DoubletFinder') Download and installation of DoubletFinder
library(Seurat)
library(ggplot2)
library(tidyverse)
library(DoubletFinder)

# create counts matrix
cts <- ReadMtx(mtx = '~/Desktop/scRNA-seq_Brush1/filtered_feature_bc_matrix/matrix.mtx.gz',
               features = '~/Desktop/scRNA-seq_Brush1/filtered_feature_bc_matrix/features.tsv.gz',
               cells = '~/Desktop/scRNA-seq_Brush1/filtered_feature_bc_matrix/barcodes.tsv.gz')

cts[1:10,1:10]

# create Seurat object
brush1.seurat <- CreateSeuratObject(counts = cts)
str(brush1.seurat)

# QC and Filtering
# explore QC

brush1.seurat$mitoPercent <- PercentageFeatureSet(brush1.seurat, pattern = '^MT-')

brush1.seurat.filtered <- subset(brush1.seurat, subset = nCount_RNA > 800 &
                                 nFeature_RNA > 500 &
                                 mitoPercent < 10)

brush1.seurat
brush1.seurat.filtered

# pre-process standard workflow
brush1.seurat.filtered <- NormalizeData(object = brush1.seurat.filtered)
brush1.seurat.filtered <- FindVariableFeatures(object = brush1.seurat.filtered)
brush1.seurat.filtered <- ScaleData(object = brush1.seurat.filtered)
brush1.seurat.filtered <- RunPCA(object = brush1.seurat.filtered)
ElbowPlot(brush1.seurat.filtered)
brush1.seurat.filtered <- FindNeighbors(object = brush1.seurat.filtered, dims = 1:20)
brush1.seurat.filtered <- FindClusters(object = brush1.seurat.filtered)
brush1.seurat.filtered <- RunUMAP(object = brush1.seurat.filtered, dims = 1:20)

DimPlot(brush1.seurat.filtered, reduction = "umap", label = T)

# ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
paramSweep <- function(seu, PCs=1:10, sct = FALSE, num.cores=1) {
  require(Seurat); require(fields); require(parallel)
  ## Set pN-pK param sweep ranges
  pK <- c(0.0005, 0.001, 0.005, seq(0.01,0.3,by=0.01))
  pN <- seq(0.05,0.3,by=0.05)
  
  ## Remove pK values with too few cells
  min.cells <- round(nrow(seu@meta.data)/(1-0.05) - nrow(seu@meta.data))
  pK.test <- round(pK*min.cells)
  pK <- pK[which(pK.test >= 1)]
  
  ## Extract pre-processing parameters from original data analysis workflow
  orig.commands <- seu@commands
  
  ## Down-sample cells to 10000 (when applicable) for computational effiency
  if (nrow(seu@meta.data) > 10000) {
    real.cells <- rownames(seu@meta.data)[sample(1:nrow(seu@meta.data), 10000, replace=FALSE)]
    data <- seu@assays$RNA$counts[ , real.cells]
    n.real.cells <- ncol(data)
  }
  
  if (nrow(seu@meta.data) <= 10000){
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA$counts
    n.real.cells <- ncol(data)
  }
  
  ## Iterate through pN, computing pANN vectors at varying pK
  #no_cores <- detectCores()-1
  if(num.cores>1){
    require(parallel)
    cl <- makeCluster(num.cores)
    output2 <- mclapply(as.list(1:length(pN)),
                        FUN = parallel_paramSweep,
                        n.real.cells,
                        real.cells,
                        pK,
                        pN,
                        data,
                        orig.commands,
                        PCs,
                        sct,mc.cores=num.cores)
    stopCluster(cl)
  }else{
    output2 <- lapply(as.list(1:length(pN)),
                      FUN = parallel_paramSweep,
                      n.real.cells,
                      real.cells,
                      pK,
                      pN,
                      data,
                      orig.commands,
                      PCs,
                      sct)
  }
  
  ## Write parallelized output into list
  sweep.res.list <- list()
  list.ind <- 0
  for(i in 1:length(output2)){
    for(j in 1:length(output2[[i]])){
      list.ind <- list.ind + 1
      sweep.res.list[[list.ind]] <- output2[[i]][[j]]
    }
  }
  
  ## Assign names to list of results
  name.vec <- NULL
  for (j in 1:length(pN)) {
    name.vec <- c(name.vec, paste("pN", pN[j], "pK", pK, sep = "_" ))
  }
  names(sweep.res.list) <- name.vec
  return(sweep.res.list)
  
}
sweep.res.list_brush1 <- paramSweep(brush1.seurat.filtered, PCs = 1:20, sct = FALSE) # This function was changed as on the lines up. Before running this line, run the function.
sweep.stats_brush1 <- summarizeSweep(sweep.res.list_brush1, GT = FALSE)
bcmvn_brush1 <- find.pK(sweep.stats_brush1)

ggplot(bcmvn_brush1, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- bcmvn_brush1 %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- brush1.seurat.filtered@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.076*nrow(brush1.seurat.filtered@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Need to run this function to make the rest of the code work
doubletFinder_v3_SeuratV5 <- function (seu, PCs, pN = 0.25, pK, nExp, reuse.pANN = FALSE, 
                                       sct = FALSE, annotations = NULL) 
{
  require(Seurat)
  require(fields)
  require(KernSmooth)
  if (reuse.pANN != FALSE) {
    pANN.old <- seu@meta.data[, reuse.pANN]
    classifications <- rep("Singlet", length(pANN.old))
    classifications[order(pANN.old, decreasing = TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("DF.classifications", pN, pK, nExp, 
                          sep = "_")] <- classifications
    return(seu)
  }
  if (reuse.pANN == FALSE) {
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA$counts[, real.cells]
    n_real.cells <- length(real.cells)
    n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)
    print(paste("Creating", n_doublets, "artificial doublets...", 
                sep = " "))
    real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
    real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
    doublets <- (data[, real.cells1] + data[, real.cells2])/2
    colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
    data_wdoublets <- cbind(data, doublets)
    if (!is.null(annotations)) {
      stopifnot(typeof(annotations) == "character")
      stopifnot(length(annotations) == length(Cells(seu)))
      stopifnot(!any(is.na(annotations)))
      annotations <- factor(annotations)
      names(annotations) <- Cells(seu)
      doublet_types1 <- annotations[real.cells1]
      doublet_types2 <- annotations[real.cells2]
    }
    orig.commands <- seu@commands
    if (sct == FALSE) {
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
      print("Normalizing Seurat object...")
      seu_wdoublets <- NormalizeData(seu_wdoublets, normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method, 
                                     scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor, 
                                     margin = orig.commands$NormalizeData.RNA@params$margin)
      print("Finding variable genes...")
      seu_wdoublets <- FindVariableFeatures(seu_wdoublets, 
                                            selection.method = orig.commands$FindVariableFeatures.RNA$selection.method, 
                                            loess.span = orig.commands$FindVariableFeatures.RNA$loess.span, 
                                            clip.max = orig.commands$FindVariableFeatures.RNA$clip.max, 
                                            mean.function = orig.commands$FindVariableFeatures.RNA$mean.function, 
                                            dispersion.function = orig.commands$FindVariableFeatures.RNA$dispersion.function, 
                                            num.bin = orig.commands$FindVariableFeatures.RNA$num.bin, 
                                            binning.method = orig.commands$FindVariableFeatures.RNA$binning.method, 
                                            nfeatures = orig.commands$FindVariableFeatures.RNA$nfeatures, 
                                            mean.cutoff = orig.commands$FindVariableFeatures.RNA$mean.cutoff, 
                                            dispersion.cutoff = orig.commands$FindVariableFeatures.RNA$dispersion.cutoff)
      print("Scaling data...")
      seu_wdoublets <- ScaleData(seu_wdoublets, features = orig.commands$ScaleData.RNA$features, 
                                 model.use = orig.commands$ScaleData.RNA$model.use, 
                                 do.scale = orig.commands$ScaleData.RNA$do.scale, 
                                 do.center = orig.commands$ScaleData.RNA$do.center, 
                                 scale.max = orig.commands$ScaleData.RNA$scale.max, 
                                 block.size = orig.commands$ScaleData.RNA$block.size, 
                                 min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block)
      print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets, features = orig.commands$ScaleData.RNA$features, 
                              npcs = length(PCs), rev.pca = orig.commands$RunPCA.RNA$rev.pca, 
                              weight.by.var = orig.commands$RunPCA.RNA$weight.by.var, 
                              verbose = FALSE)
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[, 
                                                                PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets)
      gc()
    }
    if (sct == TRUE) {
      require(sctransform)
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
      print("Running SCTransform...")
      seu_wdoublets <- SCTransform(seu_wdoublets)
      print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs))
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[, 
                                                                PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets)
      gc()
    }
    print("Calculating PC distance matrix...")
    dist.mat <- fields::rdist(pca.coord)
    print("Computing pANN...")
    pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, 
                                 ncol = 1))
    if (!is.null(annotations)) {
      neighbor_types <- as.data.frame(matrix(0L, nrow = n_real.cells, 
                                             ncol = length(levels(doublet_types1))))
    }
    rownames(pANN) <- real.cells
    colnames(pANN) <- "pANN"
    k <- round(nCells * pK)
    for (i in 1:n_real.cells) {
      neighbors <- order(dist.mat[, i])
      neighbors <- neighbors[2:(k + 1)]
      pANN$pANN[i] <- length(which(neighbors > n_real.cells))/k
      if (!is.null(annotations)) {
        for (ct in unique(annotations)) {
          neighbor_types[i, ] <- table(doublet_types1[neighbors - 
                                                        n_real.cells]) + table(doublet_types2[neighbors - 
                                                                                                n_real.cells])
          neighbor_types[i, ] <- neighbor_types[i, ]/sum(neighbor_types[i, 
          ])
        }
      }
    }
    print("Classifying doublets..")
    classifications <- rep("Singlet", n_real.cells)
    classifications[order(pANN$pANN[1:n_real.cells], decreasing = TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("pANN", pN, pK, nExp, sep = "_")] <- pANN[rownames(seu@meta.data), 
                                                                    1]
    seu@meta.data[, paste("DF.classifications", pN, pK, nExp, 
                          sep = "_")] <- classifications
    if (!is.null(annotations)) {
      colnames(neighbor_types) = levels(doublet_types1)
      for (ct in levels(doublet_types1)) {
        seu@meta.data[, paste("DF.doublet.contributors", 
                              pN, pK, nExp, ct, sep = "_")] <- neighbor_types[, 
                                                                              ct]
      }
    }
    return(seu)
  }
}
# run doubletFinder 
brush1.seurat.filtered <- doubletFinder_v3_SeuratV5(brush1.seurat.filtered, 
                                         PCs = 1:20, 
                                         pN = 0.25, 
                                         pK = pK, 
                                         nExp = nExp_poi.adj,
                                         reuse.pANN = FALSE, sct = FALSE)


# visualize doublets
DimPlot(brush1.seurat.filtered, reduction = 'umap', group.by = "DF.classifications_0.25_0.22_630") # This element DF.classifications_0.25_0.22_630 changes every time. Go find it in @meta.data


# number of singlets and doublets
table(brush1.seurat.filtered@meta.data$DF.classifications_0.25_0.21_691)







# 4. run scDblFinde to identify doublets, SingleCellExperiment object
Idents(brush1.seurat.filtered)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("scDblFinder")
install.packages("scater")
library(scater)
suppressPackageStartupMessages(library(scDblFinder))
sce <- scDblFinder(GetAssayData(brush1.seurat.filtered, slot= "counts"), clusters = Idents(brush1.seurat.filtered))










