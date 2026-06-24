# Script to integrate across conditions using Harmony and some initial analysis
# NB It starts with a seurat object created



setwd("~/Desktop/scRNA-seq_Brush1/IRIS/")

# set seed for reproducibility
#set.seed(1234)

library(harmony)
library(Seurat)
library(SeuratData)
library(tidyverse)
library(ggplot2)

# Load the data
merged_seurat <- LoadSeuratRds("~/Desktop/scRNA-seq_Brush1/IRIS/merged_seurat_obj.rds")


# Look at the Seurat object
str(merged_seurat)


# QC and filtering
merged_seurat$mito.percent <- PercentageFeatureSet(merged_seurat, pattern = '^MT-')
View(merged_seurat@meta.data)
# explore QC

# filter
merged_seurat
merged_seurat.filtered <- subset(merged_seurat, subset = nCount_RNA > 800 &
                          nFeature_RNA > 200 & 
                          mito.percent < 10) # <10 suggested by Osorio D, Cai JJ. 
                                            # Systematic determination of the mitochondrial proportion in human 
                                            # and mice tissues for single-cell RNA-sequencing data quality control. 
                                            # Bioinformatics. 2021 May 17;37(7):963-967. 
                                            # doi: 10.1093/bioinformatics/btaa751. PMID: 32840568; PMCID: PMC8599307.
merged_seurat.filtered

# standard workflow steps
merged_seurat.filtered <- NormalizeData(merged_seurat.filtered)
merged_seurat.filtered <- FindVariableFeatures(merged_seurat.filtered)
merged_seurat.filtered <- ScaleData(merged_seurat.filtered)
merged_seurat.filtered <- RunPCA(merged_seurat.filtered)
ElbowPlot(merged_seurat.filtered)
merged_seurat.filtered <- RunUMAP(merged_seurat.filtered, dims = 1:20, reduction = 'pca')

before <- DimPlot(merged_seurat.filtered, reduction = 'umap', group.by = 'Type')


# run Harmony -----------
merged_seurat.harmony <- merged_seurat.filtered %>%
  RunHarmony(group.by.vars = 'Type', plot_convergence = FALSE)

merged_seurat.harmony@reductions

merged_seurat.harmony.embed <- Embeddings(merged_seurat.harmony, "harmony")
merged_seurat.harmony.embed[1:10,1:10]

# Do UMAP and clustering using ** Harmony embeddings instead of PCA **
merged_seurat.harmony <- merged_seurat.harmony %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.5)

# visualize 
after <- DimPlot(merged_seurat.harmony, reduction = 'umap', group.by = 'Type')

before|after

# Save data
saveRDS(merged_seurat.harmony, "~/Desktop/scRNA-seq_Brush1/IRIS/merged_seurat.harmony_mt10.rds")


# Downstream analysis
# script to identify cluster identity -----------------
# Finding markers in every cluster
# Finding conserved markers 
# Finding markers DE between conditions
library(Seurat)
library(tidyverse)

# load data
merged_seurat.harmony <- readRDS('~/Desktop/scRNA-seq_Brush1/IRIS/merged_seurat.harmony_mt10.rds')
str(merged_seurat.harmony)
View(merged_seurat.harmony@meta.data)

# visualize data
clusters <- DimPlot(merged_seurat.harmony, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE, raster=FALSE)
condition <- DimPlot(merged_seurat.harmony, reduction = 'umap', group.by = 'Type')
condition|clusters

# findAll markers -----------------

FindAllMarkers(merged_seurat.harmony,      # This is to be used to find the cell type of each cluster. It compares one cluster to all the others
               logfc.threshold = 0.25,
               min.pct = 0.1,
               only.pos = TRUE,
               test.use = 'DESeq2',
               slot = 'counts')


# findConserved markers ------------- This is better when you have to compare conditions. It's a 1 to 1 comparison

# Notes:
# slot depends on the type of the test used, 
# default is data slot that stores normalized data
# DefaultAssay(merged_seurat.harmony) <- 'RNA'

DefaultAssay(merged_seurat.harmony)

library(metap)
library(BiocManager)
library(multtest)

merged_seurat.harmony <- JoinLayers(object = merged_seurat.harmony)

markers_cluster6 <- FindConservedMarkers(merged_seurat.harmony,
                                         ident.1 = "6" ,
                                         grouping.var = 'Type')

head(markers_cluster6)
write.table(markers_cluster1, "~/Desktop/scRNA-seq_Brush1/IRIS/markers_cluster6.csv", sep = ",")

# let's visualize top features
FeaturePlot(merged_seurat.harmony, features = c("CYP27A1",
                                                "MARCO",
                                                "FABP4"), min.cutoff = 'q10')

# Markers from Sikkema, L., Ramírez-Suástegui, C., Strobl, D.C. et al. An integrated cell atlas of the lung in health and disease. Nat Med 29, 1563–1577 (2023). https://doi.org/10.1038/s41591-023-02327-2
# table downloaded - 41591_2023_2327_MOESM3_ESM

DotPlot(merged_seurat.harmony, features = c("FXYD3","EPCAM","ELF3","IGFBP2","SERPINF1","TSPAN1")) + RotatedAxis() # Airway_epithelium markers

DotPlot(merged_seurat.harmony, features = c("KRT15","KRT17","TP63", "KRT5", "DLK2")) + RotatedAxis() # Basal_resting markers

DotPlot(merged_seurat.harmony, features = c("KRT5","SERPINB3")) + RotatedAxis() # Suprabasal markers

DotPlot(merged_seurat.harmony, features = c("CYP2F1","SCGB3A1","BPIFB1")) + RotatedAxis() # Club_non-nasal markers

DotPlot(merged_seurat.harmony, features = c("CYP27A1","MARCO","FABP4")) + RotatedAxis() # Alveolar_macrophages markers

DotPlot(merged_seurat.harmony, features = c("MUC5AC", "MUC5B")) + RotatedAxis() # Secretory cells (Goblet + Club) markers 

DotPlot(merged_seurat.harmony, features = c("FOXJ1", "TPPP3", "SNTN")) + RotatedAxis() # Ciliated cells markers

DotPlot(merged_seurat.harmony, features = c("BSND","IGF1","CLCNKB")) + RotatedAxis() # Ionocytes cells markers



DotPlot(merged_seurat.harmony, features = c("BSND",
                                            "IGF1",
                                            "CLCNKB")) + RotatedAxis()

# 
FeaturePlot(merged_seurat.harmony, features = c("MUC5AC", "MUC5B"),raster=FALSE) + RotatedAxis() 

# Markers seleted by me
DotPlot(merged_seurat.harmony, features = c("TP63", "KRT5", "KRT6C", "KRT14", "KRT15")) + RotatedAxis() # Basal cells markers
DotPlot(merged_seurat.harmony, features = c("KRT4", "KRT13", "KRT16", "KRT23")) + RotatedAxis() # Suprabasal cell markers
DotPlot(merged_seurat.harmony, features = c("SCGB1A1", "KRT5", "SCGB3A2", "CYP2F1", "CCDC50")) + RotatedAxis() # Club cells markers
DotPlot(merged_seurat.harmony, features = c("SPDEF", "FOXQ1", "MUC5AC", "FOXA3", "CLCA1" )) + RotatedAxis() # Goblet cells markers
DotPlot(merged_seurat.harmony, features = c("DNAH5", "SPEF2", "PIFO", "FOXJ1", "DNAI1", "CCDC40", "RSPH1", "TEKT1", "TUBB4B", "SNTN" )) + RotatedAxis() # Ciliated cells markers
DotPlot(merged_seurat.harmony, features = c("LINC01187", "ATP6V1G3", "FOXI1", "TMPRSS11E", "BSND", "SFTPB", "CFTR")) + RotatedAxis() # Ionocytes cells markers
DotPlot(merged_seurat.harmony, features = c("ZNF804A","CD86", "CD300LB", "LILRA2", "CD83", "ITGAX" )) + RotatedAxis() # Dendritic cells and Macrophages markers                                            



# Heatmaps clusters markers
clusters.markers <- FindMarkers(merged_seurat.harmony, ident.1 = 'Basal cells', logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
clusters.markers <- FindAllMarkers(merged_seurat.harmony, only.pos = TRUE)
clusters.markers%>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(merged_seurat.harmony, features = top10$gene) + NoLegend()



                                            
                                            
                                            

# min-cut off explanation:
seq(1,5)
SetQuantile('q50', seq(1,5))
SetQuantile('q10', seq(1,5))


# To calculate the number of cells per cluster from https://satijalab.org/seurat/archive/v3.0/interaction_vignette
table(Idents(merged_seurat.harmony))
# What proportion of cells are in each cluster?
prop.table(table(Idents(merged_seurat.harmony)))
# Number of cells per cluster per identity (Type)
cell_prop_per_condition <- table(merged_seurat.harmony@active.ident, merged_seurat.harmony@meta.data$Type)

cell_prop_per_condition_df <- as.data.frame(cell_prop_per_condition)
names(cell_prop_per_condition_df)[1] <- "Cell_type"
names(cell_prop_per_condition_df)[2] <- "Condition"


ggplot(cell_prop_per_condition_df, aes(fill=Cell_type, y=Freq, x=Condition)) + 
  geom_bar(position="fill", stat="identity")




# rename cluster ident
Idents(merged_seurat.harmony)
merged_seurat.harmony <- RenameIdents(merged_seurat.harmony, `13` = 'Ionocytes')


DimPlot(merged_seurat.harmony, reduction = 'umap', label = T)

# cells already have annotations provided in the metadata
View(merged_seurat.harmony@meta.data)

# Settings cluster identities is an iterative step
# multiple approaches could be taken - automatic/manual anotations (sometimes both)
# need to make sure each cell type forms a separate cluster

# setting Idents as Seurat annotations provided (also a sanity check!)
Idents(merged_seurat.harmony) <- merged_seurat.harmony@meta.data$seurat_annotations
Idents(merged_seurat.harmony)

DimPlot(merged_seurat.harmony, reduction = 'umap', label = TRUE)


# findMarkers between conditions ---------------------
merged_seurat.harmony$celltype.cnd <- paste0(merged_seurat.harmony$seurat_annotations,'_', merged_seurat.harmony$Type)
View(merged_seurat.harmony@meta.data)
Idents(merged_seurat.harmony) <- merged_seurat.harmony$celltype.cnd

DimPlot(merged_seurat.harmony, reduction = 'umap', label = TRUE)

# find markers
b.interferon.response <- FindMarkers(merged_seurat.harmony, ident.1 = 'CD16 Mono_Type', ident.2 = 'CD16 Mono_CTRL')

head(b.interferon.response)

# plotting conserved features vs DE features between conditions
head(markers_cluster3)

FeaturePlot(merged_seurat.harmony, features = c('MUC5AC'), split.by = 'Type', min.cutoff = 'q10')

# Integrating stimulated vs. control to learn cell-type specific responses (https://satijalab.org/seurat/archive/v3.1/immune_alignment.html)

DimPlot(merged_seurat.harmony, reduction = "umap", split.by = "Type")







##### - This is to compare the expression between A and B in different cell types
library(ggplot2)
library(cowplot)

theme_set(theme_cowplot())
Basal <- subset(merged_seurat.harmony, idents = "Basal cells")
Idents(Basal) <- "Type"
avg.Basal <- log1p(AverageExpression(Basal, verbose = FALSE)$RNA)
#avg.Basal$gene <- rownames(avg.Basal)
avg.Basal@Dimnames[[1]] <- rownames(avg.Basal)

p1 <- ggplot(avg.Basal, aes(A, B)) + geom_point() + ggtitle("Basal cells")
p1
###
theme_set(theme_cowplot())
t.cells <- subset(merged_seurat.harmony, idents = "T-cells")
Idents(t.cells) <- "Type"
avg.t.cells <- log1p(AverageExpression(t.cells, verbose = FALSE)$RNA)
#avg.Basal$gene <- rownames(avg.Basal)
avg.t.cells@Dimnames[[1]] <- rownames(avg.t.cells)

p2 <- ggplot(avg.t.cells, aes(A, B)) + geom_point() + ggtitle("T-cells")
p2
###
theme_set(theme_cowplot())
cil.cells <- subset(merged_seurat.harmony, idents = "Ciliated cells")
Idents(cil.cells) <- "Type"
avg.cil.cells <- log1p(AverageExpression(cil.cells, verbose = FALSE)$RNA)
#avg.Basal$gene <- rownames(avg.Basal)
avg.cil.cells@Dimnames[[1]] <- rownames(avg.cil.cells)

p3 <- ggplot(avg.cil.cells, aes(A, B)) + geom_point() + ggtitle("Ciliated cells")
p3
###


# # DESeq2 analysis
# 
# # Extract raw counts and metadata to create SingleCellExperiment object
# counts <- merged_seurat.harmony@assays$RNA@layers
# 
# metadata <- merged_seurat.harmony@meta.data
# 
# # Set up metadata as desired for aggregation and DE analysis
# metadata$cluster_id <- factor(merged_seurat.harmony@active.ident)
# 
# # Create single cell experiment object
# sce <- merged_seurat.harmony(assays = list(counts = layers), 
#                             colData = metadata)


library(DESeq2) # If you have problem with this, ensure you load DESeq2 firs. Save the workspace, close R and reopen the file and load DESeq2
library(ExperimentHub)
library(Seurat)

library(tidyverse)




# PSEUDOBULK ANALYSIS
# pseudo-bulk workflow -----------------
# Acquiring necessary metrics for aggregation across cells in a sample
# 1. counts matrix - sample level
# counts aggregate to sample level

View(merged_seurat.harmony@meta.data)

merged_seurat.harmony$samples <- paste0(merged_seurat.harmony$Type, merged_seurat.harmony$Patient)
merged_seurat.harmony$cluster_id <- paste0(factor(merged_seurat.harmony@active.ident))


DefaultAssay(merged_seurat.harmony)

cts <- AggregateExpression(merged_seurat.harmony, 
                           group.by = c("cluster_id", "samples"),
                           assays = 'RNA',
                           slot = "counts",
                           return.seurat = FALSE)

cts <- cts$RNA

# transpose
cts.t <- t(cts)


# convert to data.frame
cts.t <- as.data.frame(cts.t)

# get values where to split
splitRows <- gsub('_.*', '', rownames(cts.t))


# split data.frame
cts.split <- split.data.frame(cts.t,
                              f = factor(splitRows))

# fix colnames and transpose
cts.split.modified <- lapply(cts.split, function(x){
  rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
  t(x)
  
})

#gsub('.*_(.*)', '\\1', 'B cells_ctrl101')


############################################
## Let's run DE analysis with Basal cells ##
############################################
# 1. Get counts matrix
counts_bas.cell <- cts.split.modified$`Basal cells`

# 2. generate sample level metadata
colData <- data.frame(samples = colnames(counts_bas.cell))

colData <- colData %>%
  mutate(condition = ifelse(grepl('B', samples), 'PostMEPO', 'PreMEPO')) %>%
  column_to_rownames(var = 'samples')

# get more information from metadata




# perform DESeq2 --------
# Create DESeq2 object   
dds <- DESeqDataSetFromMatrix(countData = counts_bas.cell,
                              colData = colData,
                              design = ~ condition)

# filter
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "PreMEPO")
# run DESeq2
dds <- DESeq(dds)

summary(results(dds, alpha=0.05))

# Check the coefficients for the comparison
resultsNames(dds)


# Generate results object
res <- results(dds, name = "condition_PostMEPO_vs_PreMEPO")
res


# How many adjusted p-values were less than 0.1?
sum(res$padj < 0.1, na.rm=TRUE)

# save res file
write.csv(res, file = "~/Desktop/scRNA-seq_Brush1/IRIS/res_samp7-18_basal_corrected.csv")


# generate volcanoplot
library(EnhancedVolcano)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "DE Basal cells pre/post MEPO corrected",
                pCutoff = 0.05,
                FCcutoff = 2,
                col = c("grey", "grey", 'grey', "red3"),
                #selectLab = c("GSTA1", "GSTA2","HLA-DQA1", "HLA-DRB1", "NQO1", "B2M", "PRDX1"),
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                drawConnectors = FALSE,
                boxedLabels = FALSE,
                legendPosition = "top")

# generate MA plot
plotMA(res, ylim=c(-8,8))

# PCA
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)

pPCA <- plotPCA(vsd, intgroup=c("condition", "sizeFactor"))
pPCA

vsd@colData$type <- c("18", "7", "18", "7")

pcaData <- plotPCA(vsd, intgroup=c("condition", "type"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

# Comparison in featureplot
FeaturePlot(merged_seurat.harmony, features = c('PIGR', 'GSTA1', 'CYP2F1'), split.by = 'Type', min.cutoff = 'q9')

# Comparison by Dotplot
Idents(merged_seurat.harmony) <- factor(Idents(merged_seurat.harmony), levels = c("Basal cells"))
markers.to.plot <- c("IL5", "IL5RA")
DotPlot(merged_seurat.harmony, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8, 
        split.by = "Type") + RotatedAxis()


## Try to detect doublets with DoubletFinder 


library(DoubletFinder)


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
sweep.res.list_IRIS <- paramSweep(merged_seurat.harmony, PCs = 1:20, sct = FALSE) # This function was changed as on the lines up. Before running this line, run the function.
sweep.stats_IRIS <- summarizeSweep(sweep.res.list_IRIS, GT = FALSE)
bcmvn_IRIS <- find.pK(sweep.stats_IRIS)

ggplot(bcmvn_IRIS, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- bcmvn_IRIS %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- merged_seurat.harmony@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.076*nrow(merged_seurat.harmony@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
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
merged_seurat.harmony <- doubletFinder_v3_SeuratV5(merged_seurat.harmony, 
                                                    PCs = 1:20, 
                                                    pN = 0.25, 
                                                    pK = pK, 
                                                    nExp = nExp_poi.adj,
                                                    reuse.pANN = FALSE, sct = FALSE)


# visualize doublets
DimPlot(merged_seurat.harmony, reduction = 'umap', group.by = "DF.classifications_0.25_0.005_164") # This element: DF.classifications_0.25_0.005_164, changes every time. Go find it in @meta.data


# number of singlets and doublets
table(merged_seurat.harmony@meta.data$DF.classifications_0.25_0.005_164)







# 4. run scDblFinde to identify doublets, SingleCellExperiment object
Idents(merged_seurat.harmony)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("scDblFinder")
install.packages("scater")
library(scater)
suppressPackageStartupMessages(library(scDblFinder))
sce <- scDblFinder(GetAssayData(merged_seurat.harmony, slot= "counts"), clusters = Idents(merged_seurat.harmony))











