library(SeuratExtend)
library(pheatmap)
library(patchwork)
library(dplyr)
library(ggplot2)
library(pheatmap)
epithelial_obj <- readRDS("/cephfs/volumes/hpc_home/k2260933/7a364635-0534-4a3e-b842-27e98e9382dc/Seurat_objs/epithelial_obj.rds")
#
# Trajectory Analysis 
# Diffusion Map Calculation
# 
# Palantir uses diffusion maps for dimensionality reduction to infer trajectories. Here's how to compute and visualize them:
# Compute diffusion map
merged_seurat.harmony<-Palantir.RunDM(epithelial_obj)
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
treatments <- unique(merged_seurat.harmony$Type)

plots <- lapply(treatments, function(trt) {
  seu_sub <- subset(merged_seurat.harmony, Type == trt)
  
  DimPlot2(
    seu_sub,
    reduction = "ms",
    group.by = "annotation1",
    label = TRUE,
    repel = TRUE
  ) + 
    ggtitle(trt)
})

wrap_plots(plots)

# Pseudotime Calculation
# Pseudotime ordering assigns each cell a time point in a trajectory, indicating its progression along a developmental path:
# Calculate pseudotime with a specified start cell
# Extract expression vector for TP63
tp63_expr <- FetchData(merged_seurat.harmony, vars = "TP63")

# Get the name of the cell with the maximum expression
top_cell <- rownames(tp63_expr)[which.max(tp63_expr$TP63)]
top_cell

merged_seurat.harmony <- Palantir.Pseudotime(merged_seurat.harmony, start_cell = "P19_A_AATAGAGTCAATCTCT-1")
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
colnames(ps)[3:6] <- c("fate1", "fate2", "fate3", "fate4")
merged_seurat.harmony@meta.data[,colnames(ps)] <- ps

# Save updated seurat object
saveRDS(merged_seurat.harmony, "~/Seurat_objs/merged_seurat.harmony.pseudotime.rds")

-----------------------------------------------------------------------------------------------------------------
  # Load dataset with pseudotime
  merged_seurat.harmony <- readRDS("./Seurat_objs/merged_seurat.harmony.pseudotime.rds")
-----------------------------------------------------------------------------------------------------------------
  ps <- merged_seurat.harmony@misc$Palantir$Pseudotime
colnames(ps)[3:6] <- c("fate1", "fate2", "fate3", "fate4")
merged_seurat.harmony@meta.data[,colnames(ps)] <- ps

# Visualize pseudotime and cell fates
DimPlot2(
  merged_seurat.harmony,
  features = colnames(ps),
  reduction = "ms",
  cols = list(continuous = "A", Entropy = "D"))

# Pseudotime comparison
plot_df <- merged_seurat.harmony@meta.data

plot_df$annotation1 <- factor(plot_df$annotation1,
                              levels = names(sort(tapply(plot_df$Pseudotime, plot_df$annotation1, median))))

ggplot(plot_df, aes(x = annotation1, y = Pseudotime, fill = Type)) +
  geom_jitter(
    aes(color = Type),
    position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9),
    alpha = 0.3,
    size = 0.5,
    show.legend = FALSE
  ) + # points first, behind
  scale_fill_manual(values = c("A" =  "grey80", "B" = "steelblue"),
                    labels = c("A" = "pre", "B" = "post")) +  #"#1f77b4", "midnightblue"
  scale_color_manual(values = c("A" =  "grey80", "B" = "steelblue"),
                     labels = c("A" = "pre", "B" = "post")) +
  geom_violin(
    position = position_dodge(width = 0.9),
    trim = FALSE,
    alpha = 0.8
  ) +  # violin on top
  theme_classic() +
  xlab("Cell type") +
  ylab("Pseudotime") +
  ggtitle("Pseudotime Distribution by Cell Type and Condition") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_blank()
  )
-----------------------------------------------------------------------------------------------------------------
  # Enhanced Dot Plots (New in v1.1.0)
  -----------------------------------------------------------------------------------------------------------------
  # Create grouped features
  grouped_features <- list(
    "Dendritic" = c("CD1C", "CLEC9A", "ZBTB46", "IRF8", "BATF3", "FLT3", "IL12B"),
    "T-cells"=c("PTPRC", "CD2"),
    "Mast cells" =c("KIT", "HPGDS", "TPSB2", "TPSAB1"),
    "Macrophages"=c("CSF1R", "CD163", "CD14","APOE", "PPARG", "MSR1", "CD68"),
    "Ionocytes"=c("SFTPB", "BSND", "TMPRSS11E", "FOXI1", "ATP6V1G3", "LINC01187","CFTR", "ASCL3"),
    "Goblet"=c("SPDEF", "MUC5B", "MUC5AC"),
    "Deuterosomal" = c("FOXN4", "DEUP1","CDC20B", "IRF8"),
    "Club"=c("CCDC50", "CYP2F1", "SCGB1A1", "CC10"),
    "Ciliated"=c("SNTN", "TUBB4B", "TEKT1", "RSPH1", "CCDC40", "DNAI1", "FOXJ1", "PIFO", "SPEF2", "DNAH5"),
    "Basal"=c("KRT15", "KRT5", "TP63","KRT17"),
    "Tuft" = c("MYB", "ASCL2", "LRMP", "RG513", "GNAT3", "GNB3", "TRPM5", "DCLK1")
  )

DotPlot2(merged_seurat.harmony, features = grouped_features, group.by = "annotation1")

# Visualization Along Trajectories
# 
# Visualizing gene expression or regulon activity along calculated trajectories can provide insights into dynamic changes:
# Create smoothed gene expression curves along trajectory
GeneTrendCurve.Palantir(
  merged_seurat.harmony,
  pseudotime.data = ps,
  # features = c("TP63", "KRT5", "SCGB1A1", "MUC5B", "FOXJ1","SMAD2", "CFTR") # Some markers of differentiation from Basal to Ciliated, Mucociliated ...
  # features = c("KRT15", "KRT5", "TP63","KRT17") # Basal
  features = c("FOXI1", "CFTR", "ASCL3", "TMPRSS11E") # Ionocytes (fate4)
  # features = c("CCDC50", "CYP2F1", "SCGB1A1") # Club
  # features = c("FOXN4", "DEUP1","CDC20B","IRF8") # Deuterosomal
  #features = c("FOXJ1", "PIFO", "SPEF2", "DNAH5") # Ciliated
) 
# To force the plots in 2 rows, save the plot as trend (trend <- GeneTrendCurve.Palatir....) and add the following
# wrap_plots(trend, ncol = 2)

# Separate A and B gene expression curve
treatments <- unique(merged_seurat.harmony$Type)

plots <- lapply(treatments, function(trt) {
  # Subset the Seurat object
  seu_sub <- merged_seurat.harmony[, merged_seurat.harmony$Type == trt]
  # Subset the pseudotime data to the same cells
  ps_sub <- ps[colnames(seu_sub), , drop = FALSE]
  # Run the plot
  GeneTrendCurve.Palantir(
    merged_seurat.harmony,
    pseudotime.data = ps,
    # features = c("TP63", "KRT5", "SCGB1A1", "MUC5B", "FOXJ1","SMAD2", "CFTR") # Some markers of differentiation from Basal to Ciliated, Mucociliated ...
    features = c("KRT15", "KRT5", "TP63","KRT17") # Basal
    # features = c("FOXI1", "CFTR", "ASCL3", "TMPRSS11E") # Ionocytes (fate4)
    # features = c("CCDC50", "CYP2F1", "SCGB1A1") # Club
    # features = c("FOXN4", "DEUP1","CDC20B","IRF8") # Deuterosomal
    # features = c("MAGI3") # Ciliated
  )})
# Combine plots
wrap_plots(plots)

# or
GeneTrendHeatmap.Palantir(
  merged_seurat.harmony,
  features = c("TP63", "KRT5", "SCGB1A1", "MUC5B", "FOXJ1", "SMAD2"),
  pseudotime.data = ps,
  lineage = "fate1"
)
# Separate A and B
treatments <- unique(merged_seurat.harmony$Type)

plots <- lapply(treatments, function(trt) {
  # Subset the Seurat object
  seu_sub <- merged_seurat.harmony[, merged_seurat.harmony$Type == trt]
  # Subset the pseudotime data to the same cells
  ps_sub <- ps[colnames(seu_sub), , drop = FALSE]
  # Run the plot
  GeneTrendHeatmap.Palantir(
    seu_sub,
    pseudotime.data = ps_sub,
    features = c("TP63", "KRT5", "SCGB1A1", "MUC5B", "FOXJ1", "SMAD2"), # markers
    # features = c("DMTF1", "PLAG1", "MBD5"),
    # features = "IL5RA",
    lineage = "fate1"
  ) + ggtitle(trt)
})

wrap_plots(plots, nrow = 1, ncol = 2)
-----------------------------------------------------------------------------------------------------------------
  # Create a gene trend heatmap for different fates TOP GENES
  -----------------------------------------------------------------------------------------------------------------  
  GeneTrendHeatmap.Palantir(
    merged_seurat.harmony,
    features = VariableFeatures(merged_seurat.harmony)[1:10],
    pseudotime.data = ps,
    lineage = "fate1"
  )

# Separate A and B
treatments <- unique(merged_seurat.harmony$Type)

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


# IL-5 Signature
IL5_sign_genes <- c("ABI1","ACVR2A","ADGRA3","ANOS1","ATP6V0D1","BCL2L11","CCDC89","CCT3","CCT5","CENPP",
                    "CFL1","CLDN16","CXXC5","DENND1B","DMTF1","DUS1L","ERCC8","FUCA2","GANC","GDA","GLB1L","KRT8",
                    "MAGI3","MBD5","MBTD1","MRPS6","ORC5","PCBD2","PDCL3","PDLIM1","PGAM1","PLAG1","RAB11FIP2",
                    "ROBO1","SCD5","ST8SIA4","STAU1","STK3","TPM2","ZCCHC7")

remodeling_genes <- c("CNN3", "TLK1", "RARB", "GSK3B", "PSMD2", "SMARCAD1", "NR3C2", "NR3C1", "CHD7", "STK3",
                      "VIM", "ARHGAP21", "ABI1", "MDK", "CFL1", "GAPDH", "CDK2AP1", "PRKD1", "MAP2K4", "SMAD2",
                      "PLTP", "TIMP1")

plots <- lapply(treatments, function(trt) {
  
  seu_sub <- merged_seurat.harmony[, merged_seurat.harmony$Type == trt]
  ps_sub <- ps[colnames(seu_sub), , drop = FALSE]
  
  p <- GeneTrendHeatmap.Palantir(
    seu_sub,
    pseudotime.data = ps_sub,
    features = remodeling_genes,
    lineage = "fate1"
  ) + ggtitle(trt)
  # Ensure consistent gene order if "Feature" is an axis
  p + scale_y_discrete(limits = remodeling_genes)
})
wrap_plots(plots)




















-----------------------------------------------------------------------------------------------------------------
  # Heatmap differences
  -----------------------------------------------------------------------------------------------------------------  
  # Define the gene order
  # IL-5 Signature
  genes <- c("ABI1","ACVR2A","ADGRA3","ANOS1","ATP6V0D1","BCL2L11","CCDC89","CCT3","CCT5","CENPP",
             "CFL1","CLDN16","CXXC5","DENND1B","DMTF1","DUS1L","ERCC8","FUCA2","GANC","GDA","GLB1L","KRT8",
             "MAGI3","MBD5","MBTD1","MRPS6","ORC5","PCBD2","PDCL3","PDLIM1","PGAM1","PLAG1","RAB11FIP2",
             "ROBO1","SCD5","ST8SIA4","STAU1","STK3","TPM2","ZCCHC7")

# Remodeling
genes <- c("CNN3", "TLK1", "RARB", "GSK3B", "PSMD2", "SMARCAD1", "NR3C2", "NR3C1", "CHD7", "STK3",
           "VIM", "ARHGAP21", "ABI1", "MDK", "CFL1", "GAPDH", "CDK2AP1", "PRKD1", "MAP2K4", "SMAD2",
           "PLTP", "TIMP1")


# ============================================================
# Heatmap comparison: Condition A vs B across pseudotime
# Marker genes: Basal and Ionocytes
# ============================================================

library(Seurat)
library(pheatmap)
library(stats)
library(grid)
library(gridExtra)

# load Surat object
merged_seurat.harmony <- readRDS("./Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/merged_seurat.harmony.pseudotime.rds")

ps <- merged_seurat.harmony@misc$Palantir$Pseudotime
colnames(ps)[3:6] <- c("fate1", "fate2", "fate3", "fate4")
merged_seurat.harmony@meta.data[,colnames(ps)] <- ps

# --- Define genes and conditions -------------------------
# Marker genes Basal cells and Ionocytes
genes <- c("KRT15", "KRT5", "TP63", "KRT17",
           "FOXI1", "CFTR", "ASCL3", "TMPRSS11E")
# Ciliated cells margers
genes <- c("DNAH5", "SPEF2", "PIFO", "FOXJ1", "DNAI1", "CCDC40", "RSPH1", "TEKT1", "TUBB4B", "SNTN")
# IL-5 Signature
genes <- c("ABI1","ACVR2A","ADGRA3","ANOS1","ATP6V0D1","BCL2L11","CCDC89","CCT3","CCT5","CENPP",
           "CFL1","CLDN16","CXXC5","DENND1B","DMTF1","DUS1L","ERCC8","FUCA2","GANC","GDA","GLB1L","KRT8",
           "MAGI3","MBD5","MBTD1","MRPS6","ORC5","PCBD2","PDCL3","PDLIM1","PGAM1","PLAG1","RAB11FIP2",
           "ROBO1","SCD5","ST8SIA4","STAU1","STK3","TPM2","ZCCHC7")
# Remodeling
genes <- c("CNN3", "TLK1", "RARB", "GSK3B", "PSMD2", "SMARCAD1", "NR3C2", "NR3C1", "CHD7", "STK3",
           "VIM", "ARHGAP21", "ABI1", "MDK", "CFL1", "GAPDH", "CDK2AP1", "PRKD1", "MAP2K4", "SMAD2",
           "PLTP", "TIMP1")

conditions <- c("A", "B")

# --- 1. Run per condition -------------------------------------------

cells_A <- colnames(merged_seurat.harmony)[merged_seurat.harmony$Type == "A"]
cells_B <- colnames(merged_seurat.harmony)[merged_seurat.harmony$Type == "B"]

p_A <- GeneTrendHeatmap.Palantir(
  merged_seurat.harmony[, cells_A],
  features        = genes,
  pseudotime.data = ps[cells_A, , drop = FALSE],
  lineage         = "fate1"
)

p_B <- GeneTrendHeatmap.Palantir(
  merged_seurat.harmony[, cells_B],
  features        = genes,
  pseudotime.data = ps[cells_B, , drop = FALSE],
  lineage         = "fate1"
)

# --- 2. Extract and reshape to matrix --------------------------------

df_to_mat <- function(df) {
  mat <- dcast(df, id ~ variable, value.var = "value")
  rownames(mat) <- mat$id
  mat$id <- NULL
  as.matrix(mat)
}

mat_A <- df_to_mat(p_A$data)
mat_B <- df_to_mat(p_B$data)

# --- 3. Align genes (same rows, same order) --------------------------

shared_genes <- intersect(rownames(mat_A), rownames(mat_B))
mat_A <- mat_A[shared_genes, ]
mat_B <- mat_B[shared_genes, ]

# --- 4. Min-max normalise across both conditions → [0, 1] ------------

minmax_rows_combined <- function(mat_A, mat_B) {
  mat_combined <- cbind(mat_A, mat_B)
  mat_norm <- t(apply(mat_combined, 1, function(x) {
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) == 0) return(rep(0.5, length(x)))
    (x - rng[1]) / (rng[2] - rng[1])
  }))
  n_A <- ncol(mat_A)
  n_B <- ncol(mat_B)
  list(
    A = mat_norm[, 1:n_A,                  drop = FALSE],
    B = mat_norm[, (n_A + 1):(n_A + n_B),  drop = FALSE]
  )
}

normed        <- minmax_rows_combined(mat_A, mat_B)
mat_A_norm    <- normed$A
mat_B_norm    <- normed$B
mat_diff_norm <- mat_B_norm - mat_A_norm # Shows B in respect of A

# --- 5. Colour schemes -----------------------------------------------
# Expression: reproduce Palantir magma-like palette

color_expr  <- colorRampPalette(c("#440154FF", "#414487FF", "#2A788EFF",
                                  "#22A884FF", "#7AD151FF", "#FDE725FF"))(100)

color_diff  <- colorRampPalette(c("#2166AC", "#F7F7F7", "#D6604D"))(100)
# color_diff  <- colorRampPalette(c("#440154FF", "#414487FF", "#2A788EFF",
#                                   "#22A884FF", "#7AD151FF", "#FDE725FF"))(100)
diff_breaks <- seq(-1, 1, length.out = 101)

# --- 6. Heatmap helper -----------------------------------------------

get_pheatmap_grob <- function(mat, title, color, breaks = NA) {
  p <- pheatmap(
    mat,
    cluster_rows  = FALSE,
    cluster_cols  = FALSE,
    color         = color,
    breaks        = breaks,
    main          = title,
    fontsize_row  = 7,
    show_colnames = FALSE,
    #gaps_row      = 4, # This introduces a gap after the 4th row
    silent        = TRUE
  )
  return(p$gtable)
}

# --- 7. Plot ---------------------------------------------------------

grob_A    <- get_pheatmap_grob(mat_A_norm,    "pre",
                               color_expr)
grob_B    <- get_pheatmap_grob(mat_B_norm,    "post",
                               color_expr)
grob_diff <- get_pheatmap_grob(mat_diff_norm, "Difference (post - pre)",
                               color_diff, breaks = diff_breaks)

grid.arrange(grob_A, grob_B, grob_diff, ncol = 3)



"#FFFFFF", "#D7E8F0", "#AFD1E1", "#87BAD2", "#5FA3C3", "#4682B4"

"#000004", "#3B0F70", "#8C2981","#DE4968", "#FE9F6D", "#FCFDBF"

"#440154FF", "#414487FF", "#2A788EFF","#22A884FF", "#7AD151FF", "#FDE725FF"

"#313695", "#FFFFFF", "#A50026"

"#00897B", "#FAFAFA", "#E64A19"

"#762A83", "#F7F7F7", "#1B7837"

"#1A237E", "#ECEFF1", "#F9A825"

"#1D78B4", "#74C4E8", "#FFFFFF", "#F4A95A", "#C94B27"
