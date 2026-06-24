# SeuratExtend - https://huayc09.github.io/SeuratExtend/
# Hua, Y., Weng, L., Zhao, F., and Rambow, F. (2025). SeuratExtend: streamlining single-cell RNA-seq analysis through an integrated and intuitive framework. Gigascience 14, giaf076. https://doi.org/10.1093/gigascience/giaf076.

library(Seurat)
library(SeuratExtend)

# load data
merged_seurat.harmony <- readRDS("~/Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/merged_seurat.harmonyCORRECTED.rds")
merged_seurat.harmony <- readRDS("~/Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/epithelial_obj.rds") # EPTITHELIAL only
# Visualizing cell clusters using DimPlot2
DimPlot2(merged_seurat.harmony, label = TRUE, theme = NoLegend(), group.by = "annotation1") + NoLegend() + theme(axis.text=element_text(size=12, face = "bold"),
                                                                                                                                    axis.title=element_text(size=12,face="bold"))

DimPlot2(merged_seurat.harmony, label = TRUE, theme = NoLegend(), group.by = "annotation1", split.by = "Type") + NoLegend() + theme(axis.text=element_text(size=12, face = "bold"),
                                                                                                                 axis.title=element_text(size=12,face="bold"))
                                                                                                                                    


# Analyzing Cluster Distribution
# 
# To check the percentage of each cluster within different samples:
# Cluster distribution bar plot
ClusterDistrBar(merged_seurat.harmony$Type, merged_seurat.harmony$annotation1, flip = FALSE) + theme(axis.text=element_text(size=12, face = "bold"),
                                                                                                     axis.title=element_text(size=12,face="bold"))

# Separate the Seurat object according to the cell Origin (Epithelial or Immune, see SeuratFull code)
# Create one Seurat object with only Epithelial cells
seurat_epithelial <- subset(
  merged_seurat.harmony,
  subset = Origin == "Epithelial"
)

# Create another Seurat object with only Immune cells
seurat_immune <- subset(
  merged_seurat.harmony,
  subset = Origin == "Immune cells"
)

saveRDS(seurat_epithelial, "~/Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/epithelial_obj.rds")
saveRDS(seurat_immune, "~/Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/immune_obj.rds")

# Marker Gene Analysis with Heatmap
# 
# To examine the marker genes of each cluster and visualize them using a heatmap:
# Calculating z-scores for variable features
genes.zscore <- CalcStats(
  merged_seurat.harmony,
  features = VariableFeatures(merged_seurat.harmony),
  group.by = "annotation1",
  order = "p",
  n = 4)
#
# Displaying heatmap
Heatmap(genes.zscore, lab_fill = "zscore")

# Enhanced Dot Plots (New in v1.1.0)
# Create grouped features
grouped_features <- list(
  "Dendritic" = c("CD1C", "CLEC9A", "ZBTB46", "IRF8", "BATF3", "FLT3", "IL12B"),
  "T-cells"=c("PTPRC", "CD2"),
    "Mast cells" =c("KIT", "HPGDS", "TPSB2", "TPSAB1"),
    "Macrophages"=c("CSF1R", "CD163", "CD14","APOE", "PPARG", "MSR1", "CD68"),
    "Ionocytes"=c("SFTPB", "BSND", "TMPRSS11E", "FOXI1", "ATP6V1G3", "LINC01187","CFTR"),
    "Goblet"=c("SPDEF", "MUC5B", "MUC5AC"),
    "Deuterosomal" = c("FOXN4", "DEUP1","CDC20B", "IRF8"),
    "Club"=c("CCDC50", "CYP2F1", "SCGB1A1"),
   "Ciliated"=c("SNTN", "TUBB4B", "TEKT1", "RSPH1", "CCDC40", "DNAI1", "FOXJ1", "PIFO", "SPEF2", "DNAH5"),
    "Basal"=c("KRT15", "KRT5", "TP63")
)
#
DotPlot2(merged_seurat.harmony, features = grouped_features, group.by = "annotation1")

# Enhanced Visualization of Marker Genes
# 
# For visualizing specific markers via a violin plot that incorporates box plots, median lines, and performs statistical testing:
# Specifying genes and cells of interest
genes <- c("TP63", "CFAP43", "SERPINB3")
cells <- merged_seurat.harmony$annotation1 %in% c("Ciliated", "Basal", "Club")
#
# Violin plot with statistical analysis
VlnPlot2(
  merged_seurat.harmony,
  features = genes,
  group.by = "annotation1",
  cells = cells,
  stat.method = "wilcox.test")

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
# DimPlot2(merged_seurat.harmony, reduction = "ms", group.by = "annotation1", label = TRUE, replel = TRUE)
DimPlot2(
  merged_seurat.harmony,
  reduction = "ms",
  group.by = "annotation1",
  label = TRUE,
  repel = TRUE
) + 
  ggtitle("Diffusion map epithelial") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Separate A and B
library(patchwork)
treatments <- unique(merged_seurat.harmony$Type)

plots <- lapply(treatments, function(trt) {
  seu_sub <- subset(merged_seurat.harmony, Type == trt)
  
  DimPlot2(
    seu_sub,
    reduction = "ms",
    group.by = "annotation1",
    label = TRUE,
    repel = TRUE
  ) + ggtitle(trt)
})

wrap_plots(plots)

# Pseudotime Calculation
# Pseudotime ordering assigns each cell a time point in a trajectory, indicating its progression along a developmental path:
# Calculate pseudotime with a specified start cell

# Selecet the cell with higher expression of TP63
# Extract expression vector for TP63
tp63_expr <- FetchData(merged_seurat.harmony, vars = "TP63")

# Get the name of the cell with the maximum expression
top_cell <- rownames(tp63_expr)[which.max(tp63_expr$TP63)]
top_cell

# Use the output of the previous to select the Starting Basal cell
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
colnames(ps)[3:5] <- c("fate1", "fate2", "fate3")
merged_seurat.harmony@meta.data[,colnames(ps)] <- ps

# Visualize pseudotime and cell fates
DimPlot2(
  merged_seurat.harmony,
  features = colnames(ps),
  reduction = "ms",
  cols = list(continuous = "A", Entropy = "D"),
  theme = NoAxes())

# Pseudotime comparison
plot_df <- merged_seurat.harmony@meta.data

plot_df$annotation1 <- factor(plot_df$annotation1,
                              levels = names(sort(tapply(plot_df$Pseudotime, plot_df$annotation1, median))))

gggplot(plot_df, aes(x = annotation1, y = Pseudotime, fill = Type)) +
  geom_jitter(
    aes(color = Type),
    position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9),
    alpha = 0.5,
    size = 0.5,
    show.legend = FALSE
  ) + # points first, behind
  scale_fill_manual(values = c("A" = "skyblue", "B" = "salmon")) +
  scale_color_manual(values = c("A" = "skyblue", "B" = "salmon")) +
  geom_violin(
    position = position_dodge(width = 0.9),
    trim = FALSE,
    alpha = 0.8
  ) +  # violin on top
  theme_classic() +
  xlab("Cell type") +
  ylab("Pseudotime") +
  ggtitle("Pseudotime Distribution by Cell Type and Condition") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Visualization Along Trajectories
# 
# Visualizing gene expression or regulon activity along calculated trajectories can provide insights into dynamic changes:
# Create smoothed gene expression curves along trajectory
GeneTrendCurve.Palantir(
  merged_seurat.harmony,
  pseudotime.data = ps,
  features = c("TP63", "KRT5", "SCGB1A1", "MUC5B", "FOXJ1","SMAD2") # Some markers of differentiation from Basal to Ciliated, Mucociliated ...
)
# or
GeneTrendHeatmap.Palantir(
  merged_seurat.harmony,
  features = c("TP63", "KRT5", "SCGB1A1", "MUC5B", "FOXJ1", "SMAD2"),
  pseudotime.data = ps,
  lineage = "fate1"
)
# Separate A and B
treatments <- unique(merged_seurat.harmony$Type)

library(patchwork)
plots <- lapply(treatments, function(trt) {
  
  # Subset the Seurat object
  seu_sub <- merged_seurat.harmony[, merged_seurat.harmony$Type == trt]
  
  # Subset the pseudotime data to the same cells
  ps_sub <- ps[colnames(seu_sub), , drop = FALSE]
  
  # Run the plot
  GeneTrendHeatmap.Palantir(
    seu_sub,
    pseudotime.data = ps_sub,
    features = c("TP63", "KRT5", "SCGB1A1", "MUC5B", "FOXJ1", "SMAD2"),
    lineage = "fate1"
  ) + ggtitle(trt)
})

wrap_plots(plots)

# Same genes order
library(patchwork)

# 1. Define all treatments
treatments <- unique(merged_seurat.harmony$Type)

# 2. Define the gene order once — same for all treatments
gene_order <- c("TP63", "KRT5", "SCGB1A1", "MUC5B", "FOXJ1")

# 3. Loop with fixed order
plots <- lapply(treatments, function(trt) {
  
  # Subset Seurat object
  seu_sub <- merged_seurat.harmony[, merged_seurat.harmony$Type == trt]
  
  # Subset pseudotime data to matching cells
  ps_sub <- ps[colnames(seu_sub), , drop = FALSE]
  
  # Generate the heatmap using the same gene order
  GeneTrendHeatmap.Palantir(
    seu_sub,
    pseudotime.data = ps_sub,
    features = gene_order,
    lineage = "fate1"
  ) + ggtitle(trt)
})

# Create a gene trend heatmap for different fates
GeneTrendHeatmap.Palantir(
  merged_seurat.harmony,
  features = VariableFeatures(merged_seurat.harmony)[1:10],
  pseudotime.data = ps,
  lineage = "fate1"
)

View(merged_seurat.harmony@meta.data)

# Separate A and B
treatments <- unique(merged_seurat.harmony$Type)

library(patchwork)
plots <- lapply(treatments, function(trt) {
  
  # Subset the Seurat object
  seu_sub <- merged_seurat.harmony[, merged_seurat.harmony$Type == trt]
  
  # Subset the pseudotime data to the same cells
  ps_sub <- ps[colnames(seu_sub), , drop = FALSE]
  
  # Run the plot
  GeneTrendHeatmap.Palantir(
    seu_sub,
    pseudotime.data = ps_sub,
    features = VariableFeatures(merged_seurat.harmony)[1:10],
    lineage = "fate1"
  ) + ggtitle(trt)
})

wrap_plots(plots)

# Same genes order
library(patchwork)

# 1. Define all treatments
treatments <- unique(merged_seurat.harmony$Type)

# 2. Define the gene order once — same for all treatments
gene_order <- VariableFeatures(merged_seurat.harmony)[1:10]

# 3. Loop with fixed order
plots <- lapply(treatments, function(trt) {
  
  # Subset Seurat object
  seu_sub <- merged_seurat.harmony[, merged_seurat.harmony$Type == trt]
  
  # Subset pseudotime data to matching cells
  ps_sub <- ps[colnames(seu_sub), , drop = FALSE]
  
  # Generate the heatmap using the same gene order
  GeneTrendHeatmap.Palantir(
    seu_sub,
    pseudotime.data = ps_sub,
    features = gene_order,
    lineage = "fate1"
  ) + ggtitle(trt)
})


# Load data slip for Epithelial or Immune cells

epithelial_data <- readRDS("~/Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/epithelial_obj.rds")
DimPlot2(epithelial_data, reduction = 'umap', label = T, group.by = 'annotation1')
ClusterDistrBar(epithelial_data$Type, epithelial_data$annotation1) + theme(axis.text=element_text(size=12, face = "bold"),
                                                                                       axis.title=element_text(size=12,face="bold"))

immune_data <- readRDS("~/Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/immune_obj.rds")
DimPlot2(immune_data, reduction = 'umap', label = T, group.by = 'annotation1')
ClusterDistrBar(immune_data$Type, immune_data$annotation1) + theme(axis.text=element_text(size=12, face = "bold"),
                                                                           axis.title=element_text(size=12,face="bold"))




# Cell type frequency as boxplot - Epithelial cells
library(ggplot2)

# Replace 'seurat_obj' with the name of your Seurat object
# Ensure 'seurat_obj' has metadata columns 'patient', 'condition', and 'cell_type'

# Extract the metadata and active identities into a data frame
cell_data <- FetchData(epithelial_data, vars = c("Patient","Type", "annotation1"))

# Calculate the number of cells per patient, condition, and cell type
cell_counts <- cell_data %>%
  group_by(Patient, Type, annotation1) %>%
  summarise(num_cells = n(), .groups = "drop")

# Calculate the total number of cells per patient and condition             #
total_counts <- cell_counts %>%                                             #
  group_by(Patient, Type) %>%                                               #
  summarise(total_cells = sum(num_cells), .groups = "drop") 

# Plot cell counts per patient per type
total_counts$Type <- factor(total_counts$Type, levels = c("B", "A"))
total_counts$Patient <- factor(total_counts$Patient, levels = c("P2", "P5", "P72", "P8", "P9", "P11", "P13", "P15", "P16", "P18", "P19", "P20"))

ggplot(total_counts, aes(x = total_cells, y = Patient, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +   # Bar plot with dodge to compare "A" vs "B"
  geom_text(aes(label = total_cells),                 # Add text labels
            position = position_dodge(width = 0.9), # Align text with bars
            hjust = 0,                              # Adjust vertical position of labels
            size = 3) +        
  labs(
    y = "Patient",
    x = "Cell Count",
    title = "Comparison of Cell Count per Patient per Type") +
  scale_fill_manual(values = c("A" = "azure4", "B" = "palegreen3")) + # Blue shades for "A" and "B"
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"), # Center and bold title
    axis.text.x = element_text(angle = 0, hjust = 1)      # Rotate x-axis labels
  )


# Join the total counts back to the cell counts and calculate fractions     #
cell_fractions <- cell_counts %>%                                           #
  left_join(total_counts, by = c("Patient", "Type")) %>%                    #
  mutate(fraction_cells = num_cells / total_cells) %>%                      #
  select(Patient, Type, annotation1, fraction_cells)                        #

# View the result
print(cell_counts)
print(cell_fractions)
# Generate the box plot
library(ggplot2)
library(ggpubr)

# Generate box plots with shared y-axis across all cell types
ggplot(cell_fractions, aes(x = annotation1, y = fraction_cells, colour = Type)) +
  geom_boxplot(
    width = 0.5,                             # Box width
    outlier.shape = NA,
    size = 0.8,
    position = position_dodge(width = 0.7)   # Space between boxes
  ) +
  geom_jitter(
    aes(color = Type),
    alpha = 0.3,  # Transparency
    position = position_dodge(width = 0.7)   # Same dodge width to align points with boxplots
  ) +
  labs(
    x = "",  # Updated x-axis label
    y = "Fraction of Cells",  # Updated y-axis label
    title = "Fraction of Epithelial Cells per Cell Type Across Conditions"  # Updated title
  ) +
  scale_fill_manual(values = c("A" = "azure4", "B" = "red")) +  # Custom fill colors
  scale_color_manual(values = c("A" = "azure4", "B" = "black")) +  # Custom point colors
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 0.5),
    strip.text = element_text(angle = 0, hjust = 1, face = "bold"),
    legend.position = "right"
  ) +
  stat_compare_means(
    aes(group = Type),                  # Compare based on Type within each annotation1 group
    method = "wilcox.test",             # Perform Wilcoxon test
    paired = F,
    label = "p.signif",                 # Show significance stars
    label.y = max(cell_fractions$fraction_cells, na.rm = TRUE) * 1.1  # Position p-values slightly above the highest point
  )


# Cell type frequency as boxplot - Immune cells
library(ggplot2)

# Replace 'seurat_obj' with the name of your Seurat object
# Ensure 'seurat_obj' has metadata columns 'patient', 'condition', and 'cell_type'

# Extract the metadata and active identities into a data frame
cell_data <- FetchData(immune_data, vars = c("Patient","Type", "annotation1"))

# Calculate the number of cells per patient, condition, and cell type
cell_counts <- cell_data %>%
  group_by(Patient, Type, annotation1) %>%
  summarise(num_cells = n(), .groups = "drop")

# Calculate the total number of cells per patient and condition             #
total_counts <- cell_counts %>%                                             #
  group_by(Patient, Type) %>%                                               #
  summarise(total_cells = sum(num_cells), .groups = "drop") 

# Plot cell counts per patient per type
total_counts$Type <- factor(total_counts$Type, levels = c("B", "A"))
total_counts$Patient <- factor(total_counts$Patient, levels = c("P2", "P5", "P72", "P8", "P9", "P11", "P13", "P15", "P16", "P18", "P19", "P20"))

ggplot(total_counts, aes(x = total_cells, y = Patient, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +   # Bar plot with dodge to compare "A" vs "B"
  geom_text(aes(label = total_cells),                 # Add text labels
            position = position_dodge(width = 0.9), # Align text with bars
            hjust = 0,                              # Adjust vertical position of labels
            size = 3) +        
  labs(
    y = "Patient",
    x = "Cell Count",
    title = "Comparison of Cell Count per Patient per Type") +
  scale_fill_manual(values = c("A" = "azure4", "B" = "palegreen3")) + # Blue shades for "A" and "B"
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"), # Center and bold title
    axis.text.x = element_text(angle = 0, hjust = 1)      # Rotate x-axis labels
  )


# Join the total counts back to the cell counts and calculate fractions     #
cell_fractions <- cell_counts %>%                                           #
  left_join(total_counts, by = c("Patient", "Type")) %>%                    #
  mutate(fraction_cells = num_cells / total_cells) %>%                      #
  select(Patient, Type, annotation1, fraction_cells)                        #

# View the result
print(cell_counts)
print(cell_fractions)
# Generate the box plot
library(ggplot2)
library(ggpubr)

# Generate box plots with shared y-axis across all cell types
ggplot(cell_fractions, aes(x = annotation1, y = fraction_cells, colour = Type)) +
  geom_boxplot(
    width = 0.5,                             # Box width
    outlier.shape = NA,
    size = 0.8,
    position = position_dodge(width = 0.7)   # Space between boxes
  ) +
  geom_jitter(
    aes(color = Type),
    alpha = 0.3,  # Transparency
    position = position_dodge(width = 0.7)   # Same dodge width to align points with boxplots
  ) +
  labs(
    x = "",  # Updated x-axis label
    y = "Fraction of Cells",  # Updated y-axis label
    title = "Fraction of Epithelial Cells per Cell Type Across Conditions"  # Updated title
  ) +
  scale_fill_manual(values = c("A" = "azure4", "B" = "red")) +  # Custom fill colors
  scale_color_manual(values = c("A" = "grey", "B" = "grey40")) +  # Custom point colors
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 0.5),
    strip.text = element_text(angle = 0, hjust = 1, face = "bold"),
    legend.position = "right"
  ) +
  stat_compare_means(
    aes(group = Type),                  # Compare based on Type within each annotation1 group
    method = "wilcox.test",             # Perform Wilcoxon test
    paired = F,
    label = "p.signif",                 # Show significance stars
    label.y = max(cell_fractions$fraction_cells, na.rm = TRUE) * 1.1  # Position p-values slightly above the highest point
  )
