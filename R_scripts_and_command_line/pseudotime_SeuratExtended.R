# Pseudotime analysis 
# SeuratExtend - https://huayc09.github.io/SeuratExtend/
# Hua, Y., Weng, L., Zhao, F., and Rambow, F. (2025). SeuratExtend: streamlining single-cell RNA-seq analysis through an integrated and intuitive framework. Gigascience 14, giaf076. https://doi.org/10.1093/gigascience/giaf076.

library(Seurat)
library(SeuratExtend)

# load data
merged_seurat.harmony <- readRDS("~/Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/merged_seurat.harmonyCORRECTED.rds")

# Visualizing Multiple Markers on UMAP
# 
# Displaying three markers on a single UMAP, using RYB coloring for each marker:
FeaturePlot3(merged_seurat.harmony, feature.1 = "IL5RA", feature.2 = "CSF2RB", pt.size = 1, color = "ryb", order = T)

# Trajectory Analysis 
# Diffusion Map Calculation
# 
# Palantir uses diffusion maps for dimensionality reduction to infer trajectories. Here's how to compute and visualize them:
# Compute diffusion map
merged_seurat.harmony<-Palantir.RunDM(merged_seurat.harmony)
## Determing nearest neighbor graph...
# Visualize the first two diffusion map dimensions
DimPlot2(merged_seurat.harmony, reduction = "ms")

# Pseudotime Calculation
#
# Pseudotime ordering assigns each cell a time point in a trajectory, indicating its progression along a developmental path:
# Calculate pseudotime with a specified start cell
merged_seurat.harmony <- Palantir.Pseudotime(merged_seurat.harmony, start_cell = "P11_A_AACGTCATCTGGAAGG-1")
## Sampling and flocking waypoints...
## Time for determining waypoints: 0.00112607479095459 minutes
## Determining pseudotime...
## Shortest path distances using 30-nearest neighbor graph...
## Time for shortest paths: 0.014574062824249268 minutes
## Iteratively refining the pseudotime...
## Correlation at iteration 1: 1.0000
## Entropy and branch probabilities...
## Markov chain construction...
## Identification of terminal states...
## Computing fundamental matrix and absorption probabilities...
## Project results to all cells...
# Store pseudotime results in meta.data for easy plotting
ps <- merged_seurat.harmony@misc$Palantir$Pseudotime
colnames(ps)[3:4] <- c("fate1", "fate2")
merged_seurat.harmony@meta.data[,colnames(ps)] <- ps

# Visualize pseudotime and cell fates
DimPlot2(
  merged_seurat.harmony,
  features = colnames(ps),
  reduction = "ms",
  cols = list(continuous = "A", Entropy = "D"),
  theme = NoAxes())

# Visualization Along Trajectories
# 
# Visualizing gene expression or regulon activity along calculated trajectories can provide insights into dynamic changes:
# Create smoothed gene expression curves along trajectory
GeneTrendCurve.Palantir(
  merged_seurat.harmony,
  pseudotime.data = ps,
  features = c("IL5", "TGF-b")
)

# Create a gene trend heatmap for different fates
GeneTrendHeatmap.Palantir(
  merged_seurat.harmony,
  features = VariableFeatures(merged_seurat.harmony)[1:10],
  pseudotime.data = ps,
  lineage = "fate1"
)