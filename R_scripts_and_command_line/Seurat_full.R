# IRIS BB pre-post
# script to integrate scRNA-Seq datasets to correct for batch effects
setwd("~/Desktop/scRNA-seq_Brush1/IRIS/Filtered_feature_bc_matrix/")


# load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(harmony)
library(SeuratData)

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

merged_seurat <- merge(P11_A, y = c(P11_B, P13_A, P13_B, P15_A, P15_B, P16_A, P16_B, P18_A, P18_B, P19_A, P19_B, P2_A, P2_B, P20_A, P20_B,
                                     P5_A, P5_B,P72_A, P72_B, P8_A, P8_B, P9_A, P9_B),
                       add.cell.ids = ls()[4:27], #this is selecting the elements output from ls() that is each seurat object created above
                       project = 'IRIS_Br')
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

# Personalised
nFeature <- VlnPlot(merged_seurat, features = c("nFeature_RNA"), pt.size = 0, cols = "grey90", annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=5000, alpha=0.2, fill="red")) & 
  #geom_hline(yintercept = 5000, color = "red") & 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 200, ymax = 5000, fill = "palegreen", alpha = 0.2) &
  theme(legend.position = 'none')
nCount <- VlnPlot(merged_seurat, features = c("nCount_RNA"), pt.size = 0, cols = "grey90") & 
  #geom_hline(yintercept = 800, color = "red") & 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 800, ymax = Inf, fill = "palegreen", alpha = 0.2) & 
  theme(legend.position = 'none')
percentMT <- VlnPlot(merged_seurat, features = c("percent.mt"), pt.size = 0, cols = "grey90") & 
  #geom_hline(yintercept = 5, color = "red") & 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 5, fill = "palegreen", alpha = 0.2) &
  theme(legend.position = 'none')

nFeature | nCount | percentMT



# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = "grey70") & theme(legend.position = 'none')
plot2 <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = "grey70") & theme(legend.position = 'none')
plot1 + plot2

# Save the object
saveRDS(merged_seurat, "~/Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/merged_seurat_obj_full_dataset_noP3.rds")

# QC and filtering
merged_seurat$mito.percent <- PercentageFeatureSet(merged_seurat, pattern = '^MT-')
View(merged_seurat@meta.data)
# explore QC

# filter
merged_seurat <- readRDS("~/Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/merged_seurat_obj_full_dataset_noP3.rds")
merged_seurat.filtered <- subset(merged_seurat, subset = nCount_RNA > 800 &
                                   nFeature_RNA > 200 & 
                                   nFeature_RNA < 5000 & # Set an upper limit of genes expressed per cell as it can limit the inclusion of doublets or more-blets
                                   mito.percent < 5) # <10 suggested by Osorio D, Cai JJ. 
saveRDS(merged_seurat.filtered, "~/Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/merged_seurat_obj_filtered200-5000_full_dataset_noP3.rds")
# Systematic determination of the mitochondrial proportion in human 
# and mice tissues for single-cell RNA-sequencing data quality control. 
# Bioinformatics. 2021 May 17;37(7):963-967. 
# doi: 10.1093/bioinformatics/btaa751. PMID: 32840568; PMCID: PMC8599307.
merged_seurat.filtered
merged_seurat.filtered <- readRDS("~/Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/merged_seurat_obj_filtered200-5000_full_dataset_noP3.rds")
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
  RunHarmony(group.by.vars = 'Patient', plot_convergence = FALSE)

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
saveRDS(merged_seurat.harmony, "~/Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/merged_seurat.harmony200-5000_full-dataset_noP3.rds")

# load data
merged_seurat.harmony <- readRDS("~/Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/merged_seurat.harmony200-5000_full-dataset_noP3.rds")
str(merged_seurat.harmony)
View(merged_seurat.harmony@meta.data)

# visualize data
clusters <- DimPlot(merged_seurat.harmony, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE, raster=FALSE)
condition <- DimPlot(merged_seurat.harmony, reduction = 'umap', group.by = 'Type')
condition|clusters

# JoinLayers
merged_seurat.harmony <- JoinLayers(object = merged_seurat.harmony)

# findAll markers -----------------
# FindAllMarkers(merged_seurat.harmony,      # This is to be used to find the cell type of each cluster. It compares one cluster to all the others
#                logfc.threshold = 0.25,
#                min.pct = 0.1,
#                only.pos = TRUE,
#                test.use = 'DESeq2',
#                slot = 'counts')
# # Notes:
# # slot depends on the type of the test used, 
# # default is data slot that stores normalized data
# # DefaultAssay(merged_seurat.harmony) <- 'RNA'
# DefaultAssay(merged_seurat.harmony)
# 
library(metap)
library(BiocManager)
library(multtest)

merged_seurat.harmony <- JoinLayers(object = merged_seurat.harmony)

markers_cluster5 <- FindConservedMarkers(merged_seurat.harmony,
                                         ident.1 = "5" ,
                                         grouping.var = 'Type')

head(markers_cluster5)
# write.table(markers_cluster1, "~/Desktop/scRNA-seq_Brush1/IRIS/markers_cluster6.csv", sep = ",")

# Markers from Sikkema, L., Ramírez-Suástegui, C., Strobl, D.C. et al. An integrated cell atlas of the lung in health and disease. Nat Med 29, 1563–1577 (2023). https://doi.org/10.1038/s41591-023-02327-2
# table downloaded - 41591_2023_2327_MOESM3_ESM
DotPlot(merged_seurat.harmony, features = c("FXYD3","EPCAM","ELF3","IGFBP2","SERPINF1","TSPAN1")) + RotatedAxis() # Airway_epithelium markers
DotPlot(merged_seurat.harmony, features = c("CYP2F1","SCGB3A1","BPIFB1")) + RotatedAxis() # Club_non-nasal markers

# Markers seleted by me
DotPlot(merged_seurat.harmony, features = c("TP63", "KRT5", "KRT15")) + RotatedAxis() # Basal cells markers
DotPlot(merged_seurat.harmony, features = c("KRT15","KRT17","TP63", "KRT5", "DLK2")) + RotatedAxis() # Basal_resting markers
DotPlot(merged_seurat.harmony, features = c("KRT4", "KRT13", "KRT16", "KRT23")) + RotatedAxis() # Suprabasal cell markers
DotPlot(merged_seurat.harmony, features = c("MUC5AC", "MUC5B")) + RotatedAxis() # Secretory cells (Goblet + Club) markers 
DotPlot(merged_seurat.harmony, features = c("SCGB1A1", "KRT5", "SCGB3A2", "CYP2F1", "CCDC50")) + RotatedAxis() # Club cells markers
DotPlot(merged_seurat.harmony, features = c("SPDEF", "FOXQ1", "MUC5AC", "FOXA3", "CLCA1" )) + RotatedAxis() # Goblet cells markers
DotPlot(merged_seurat.harmony, features = c("DNAH5", "SPEF2", "PIFO")) + RotatedAxis() # Ciliated cells markers 1
DotPlot(merged_seurat.harmony, features = c("DNAH5", "SPEF2", "PIFO", "FOXJ1", "DNAI1", "CCDC40", "RSPH1", "TEKT1", "TUBB4B", "SNTN")) + RotatedAxis() # Ciliated cells markers
DotPlot(merged_seurat.harmony, features = c("LINC01187", "ATP6V1G3", "FOXI1", "TMPRSS11E", "BSND", "SFTPB", "CFTR")) + RotatedAxis() # Ionocytes cells markers
DotPlot(merged_seurat.harmony, features = c("TPSAB1", "TPSB2", "HPGDS", "KIT")) + RotatedAxis() # Mast cell markers
DotPlot(merged_seurat.harmony, features = c("CD2", "PTPRC")) + RotatedAxis() # T-cell markers
DotPlot(merged_seurat.harmony, features = c("ZNF804A","CD86", "CD300LB", "LILRA2", "CD83", "ITGAX" )) + RotatedAxis() # Dendritic cells and Macrophages markers
DotPlot(merged_seurat.harmony, features = c("CD1C", "CLEC9A", "ZBTB46", "IRF8", "BATF3", "FLT3", "IL12B")) + RotatedAxis() # Dendritic cells markers
DotPlot(merged_seurat.harmony, features = c("CD68", "CD14", "CSF1R", "MERTK", "MSR1", "CD163", "PPARG", "APOE")) + RotatedAxis() # Macrophages
DotPlot(merged_seurat.harmony, features = c("DEUP1", "FOXN4", "CDC20B")) + RotatedAxis() # Deuterosomal
# 
FeaturePlot(merged_seurat.harmony, features = c("SCGB1A1", "KRT5", "SCGB3A2", "CYP2F1", "CCDC50"), raster = FALSE, order = TRUE) + RotatedAxis() 

# Add name to clusters
# Extract the metadata from the Seurat object
metadata <- merged_seurat.harmony@meta.data
# Define a function or use conditional logic to label entries based on another column.
# For example, we can assign labels based on some numerical value.
metadata <- metadata %>%
  mutate(annotation1 = case_when(
    seurat_clusters %in% c(4) ~ "Basal",
    seurat_clusters %in% c(0, 2, 3, 9, 11) ~ "Ciliated",
    seurat_clusters == 12 ~ "Deuterosomal",
    seurat_clusters == 17 ~ "Ionocytes",
    seurat_clusters %in% c(13) ~ "Dendritic",
    seurat_clusters %in% c(7) ~ "Macrophages",
    #seurat_clusters %in% c(19) ~ "Alveolar_Macrophages",
    seurat_clusters %in% c(8) ~ "Goblet",
    seurat_clusters %in% c(5, 14) ~ "Club",
    seurat_clusters %in% c(1, 10, 16) ~ "T-cells",
    seurat_clusters == 15 ~ "Mast cells",
    seurat_clusters == 6 ~ "Mucociliated",
    TRUE ~ "Other"  # Default label if none of the conditions above are met
  ))

# Assign the modified metadata back to the Seurat object
merged_seurat.harmony@meta.data <- metadata

# Add cell description epithelial or immune
metadata <- merged_seurat.harmony@meta.data

metadata <- metadata %>%
  mutate(Origin = case_when(
    annotation1 %in% c("Basal", "Ciliated", "Deuterosomal", "Ionocytes", "Goblet", "Club", "Mucociliated") ~ "Epithelial",
    annotation1 %in% c("Dendritic", "Macrophages", "T-cells", "Mast cells") ~ "Immune cells",
    TRUE ~ "Other" 
  ))
merged_seurat.harmony@meta.data <- metadata


# Plot the results
DimPlot(merged_seurat.harmony, reduction = 'umap', group.by = 'annotation1', label = TRUE, raster=FALSE)

marker_genes <- c("TP63", "KRT5", "KRT15", 
                  "DNAH5", "SPEF2", "PIFO", "FOXJ1", "DNAI1", "CCDC40", "RSPH1", "TEKT1", "TUBB4B", "SNTN",
                  "SCGB1A1", "CYP2F1", "CCDC50",
                  "IRF8", "DEUP1", "FOXN4", "CDC20B",
                  "MUC5AC","MUC5B", "SPDEF",
                  "LINC01187", "ATP6V1G3", "FOXI1", "TMPRSS11E", "BSND", "SFTPB", "CFTR",
                  "CD68", "CD14", "CSF1R", "MSR1", "CD163", "PPARG", "APOE",
                  "TPSAB1", "TPSB2", "HPGDS", "KIT",
                  "CD2", "PTPRC")


# Heatmap with marker genes
DoHeatmap(merged_seurat.harmony, features = c("TP63", "KRT5", "KRT15", 
                                                "DNAH5", "SPEF2", "PIFO", "FOXJ1", "DNAI1", "CCDC40", "RSPH1", "TEKT1", "TUBB4B", "SNTN",
                                                "SCGB1A1", "CYP2F1", "CCDC50",
                                                "IRF8", "DEUP1", "FOXN4", "CDC20B",
                                                "MUC5AC","MUC5B", "SPDEF",
                                                "LINC01187", "ATP6V1G3", "FOXI1", "TMPRSS11E", "BSND", "SFTPB", "CFTR",
                                                "CD68", "CD14", "CSF1R", "MSR1", "CD163", "PPARG", "APOE",
                                                "TPSAB1", "TPSB2", "HPGDS", "KIT",
                                                "CD2", "PTPRC"), 
                                              group.by = 'annotation1', slot = "scale.data", size = 3) + scale_fill_gradientn(colors = c("aliceblue", "maroon"), na.value = "white")

DotPlot(merged_seurat.harmony, c("PTPRC", "CD2",
                                 "KIT", "HPGDS", "TPSB2", "TPSAB1",
                                 "CSF1R", "CD163", "CD14","APOE", "PPARG", "MSR1", "CD68",
                                 "SFTPB", "BSND", "TMPRSS11E", "FOXI1", "ATP6V1G3", "LINC01187","CFTR",
                                 "SPDEF", "MUC5B", "MUC5AC",
                                 "FOXN4", "DEUP1","CDC20B", "IRF8",
                                 "CCDC50", "CYP2F1", "SCGB1A1",
                                 "SNTN", "TUBB4B", "TEKT1", "RSPH1", "CCDC40", "DNAI1", "FOXJ1", "PIFO", "SPEF2", "DNAH5",
                                 "KRT15", "KRT5", "TP63"), group.by = "annotation1") + 
                                coord_flip() +
                                theme(axis.text.x = element_text(angle = 45, hjust = 1),
                                axis.title = element_blank()) +
                                scale_y_discrete(limits = rev(levels(marker_genes)))
  


# Add patient demographics
# Extract the metadata from the Seurat object
metadata <- merged_seurat.harmony@meta.data
# Define a function or use conditional logic to label entries based on another column.
# For example, we can assign labels based on some numerical value.
metadata <- metadata %>%
  mutate(sex = case_when(
    Patient %in% c("P5", "P72", "P8", "P11", "P18", "P19", "P20") ~ "M",
    Patient %in% c("P2", "P9", "P13", "P15", "P16") ~ "F",
    TRUE ~ "NA"  # Default label if none of the conditions above are met
  ))

# Assign the modified metadata back to the Seurat object
merged_seurat.harmony@meta.data <- metadata

# Save data after annotation the clusters
saveRDS(merged_seurat.harmony, "~/Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/merged_seurat.harmonyCORRECTED.rds")

# load data
merged_seurat.harmony <- readRDS("~/Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/merged_seurat.harmonyCORRECTED.rds")

# Extract the metadata from the Seurat object
metadata <- merged_seurat.harmony@meta.data
# Barplot of the number of cells per condition
metadata %>% 
  ggplot(aes(x = Type, fill = Type)) + 
  
  geom_bar() +
  
  geom_text(
    stat = "count",
    aes(label = after_stat(count)),
    vjust = -0.5,   # <-- center vertically in bar
    size = 7
  ) +
  
  scale_x_discrete(
    labels = c("A" = "Pre", "B" = "Post")
  ) +
  
  scale_fill_manual(
    values = c("A" = "grey80", "B" = "steelblue"),
    labels = c("A" = "pre", "B" = "post")
  ) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.15))  # <-- more space on top
  ) +
  
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 15)
  ) +
  
  ggtitle("Number of Cells")

# Plot nCount_RNA per condition
metadata %>% 
  ggplot(aes(color = Type, x = nCount_RNA, colour = Type)) +
  geom_density(size = 1.2) +  # Adjust thickness using size argument
  scale_color_manual(values = c("grey80", "steelblue"),
                     labels = c("A" = "pre", "B" = "post")) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
  ) +
  scale_x_log10() +
  ylab("Cell density") +
  xlab("nCount_RNA (log10)") +
  ggtitle("N RNA counts")


# Cell type frequency as boxplot
library(ggplot2)

# Ensure 'seurat_obj' has metadata columns 'patient', 'condition', and 'cell_type'

# Extract the metadata and active identities into a data frame
cell_data <- FetchData(merged_seurat.harmony, vars = c("Patient","Type", "annotation1"))

# Calculate the number of cells per patient, condition, and cell type
cell_counts <- cell_data %>%
  group_by(Patient, Type, annotation1) %>%
  summarise(num_cells = n(), .groups = "drop")


# Some patients might be missing cells of a specific cell type in any of the conditions. For this you want to know the elements listed in the cell_counts = Patients x Treatments x Cell types. 
# If anithing's missing, find what:
# test_table <- cell_counts %>% 
#   group_by(Patient) %>% 
#   summarise(annotation1) %>% 
#   count(num_cells)
# Then add to the table. For this experiment is:
missing_rows_mast <- data.frame(Patient = c("P13", "P5", "P72", "P8"),
                                Type = c("A", "B", "B", "B"),
                                annotation1 = c("Mast cells","Mast cells","Mast cells","Mast cells"),
                                num_cells = c( 0,0,0,0))
missing_rows_iono_mucocil <- data.frame(Patient = c("P18", "P72", "P18", "P18"),
                                Type = c("A", "B", "A", "B"),
                                annotation1 = c("Ionocytes", "Mucociliated", "Mucociliated", "Mucociliated" ),
                                num_cells = c(0,0,0,0))
cell_counts <- rbind(cell_counts, missing_rows_iono_mucocil, missing_rows_mast)


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
                hjust = -0.2,                              # Adjust vertical position of labels
                size = 4) +        
    labs(
      y = "Patient",
      x = "Cell Count",
      title = "Comparison of Cell Count per Patient per Type") +
    scale_fill_manual(values = c("A" = "grey80", "B" = "steelblue"),
                      labels = c("A" = "pre", "B" = "post")) + # Blue shades for "A" and "B"
    scale_x_continuous(limits = c(0, 11500)) +
    theme_minimal() +
    theme(
plot.title = element_text(hjust = 0.5, face = "bold"), # Center and bold title
panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
axis.text.x = element_text(angle = 0, hjust = 1, size = 15), # Rotate x-axis labels
axis.title.x = element_text(size = 15),
axis.text.y = element_text(size = 12),
axis.title.y = element_text(size = 15)
)

# Join the total counts back to the cell counts and calculate fractions     #       CALCULATE CELL FRACTIONS
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
ggplot(cell_fractions, aes(x = annotation1, y = fraction_cells, colour = Type)) +  # To obtain the plot with the cell number change cell_fraction with cell_counts and fraction_cells with num_cells
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
    x = "",
    y = "Fraction of Cells",                                     # Updated y-axis label in case you change to cell number
    title = "Fraction of Cells per Cell Type Across Conditions"  # Updated title in case you change to cell number
  ) +
  scale_fill_manual(values = c("A" = "grey80", "B" = "steelblue"), # Custom fill colors
                    labels = c("A" = "pre", "B" = "post"),
                    name = NULL) +  
  scale_color_manual(values = c("A" = "grey80", "B" = "steelblue"), # Custom point colors
                     labels = c("A" = "pre", "B" = "post"),
                     name = NULL) +  
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 0.9, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 15),
    strip.text = element_text(angle = 0, hjust = 1, face = "bold"),
    legend.position = "right",
  
  ) +
  stat_compare_means(
    aes(group = Type),                  # Compare based on Type within each annotation1 group
    method = "wilcox.test",             # Perform Wilcoxon test
    paired = T,
    label = "p.signif",                 # Show significance stars
    hide.ns = TRUE,
    label.y = max(cell_fractions$fraction_cells, na.rm = TRUE) * 0.5,  # Position p-values slightly above the highest point
    size = 6
  )

# Other plot style - it separates every cell type
library(ggpubr)
library(dplyr)
 
   # Compute max y per group (annotation1) to position p-value labels
   label_positions <- cell_fractions %>%
       group_by(annotation1) %>%
       summarise(label_y = max(fraction_cells, na.rm = TRUE) * 1.5)

   # Merge back to cell_fractions for use in plotting
   cell_fractions_labeled <- cell_fractions %>%
       left_join(label_positions, by = "annotation1")

   ggboxplot(
     cell_fractions_labeled,
     x = "Type",
     y = "fraction_cells",
     color = "Type",
     add = "jitter",
     facet.by = "annotation1",
     palette = c("grey80", "steelblue"),
     scales = "free_y"
   ) +
     xlab("") +
     ylab("Cell fraction") +
     
     scale_color_manual(
       values = c("A" = "grey80", "B" = "steelblue"),
       labels = c("A" = "pre", "B" = "post"),
       name = NULL
     ) +
     scale_x_discrete(
       labels = c("A" = "pre", "B" = "post")
     ) +
     scale_y_continuous(expand = expansion(mult = c(0.05, 0.30))) +
     stat_compare_means(
       aes(group = Type),
       method = "wilcox.test",
       paired = TRUE,
       label = "p.signif",
       label.x = 1.45,
       hide.ns = FALSE
     ) +
     theme(
       legend.position = c(0.85, 0.15),   # bottom-right INSIDE plot
       legend.background = element_rect(fill = "white", color = "white")
     )
   
 unique_labels <- unique(cell_fractions$annotation1)

 # Plot connecting patient pre and post
 # Load libraries
 library(tidyverse)
 library(ggpubr)
 
 # Ensure correct factor ordering
 cell_fractions$Type <- factor(cell_fractions$Type, levels = c("A", "B"))
 
 # Paired Wilcoxon test per annotation1
 stats <- cell_fractions %>%
   group_by(annotation1) %>%
   summarise(
     test_result = list(
       wilcox.test(
         fraction_cells[Type == "A"],
         fraction_cells[Type == "B"],
         paired = TRUE,
         exact = FALSE
       )
     ),
     .groups = "drop"
   ) %>%
   mutate(
     p.value = map_dbl(test_result, ~ .x$p.value)
   ) %>%
   select(annotation1, p.value)
 
 print(stats)
 
 # Merge p-values back to original data for annotation in plots
 cell_fractions <- cell_fractions %>%
   left_join(stats, by = "annotation1")
 
 # Plot: paired lines for each patient
 ggplot(cell_fractions, aes(x = Type, y = fraction_cells, group = Patient, color = Patient)) +
   geom_point(size = 1) +
   geom_line(alpha = 0.5) +
   facet_wrap(~ annotation1, scales = "free_y") +
   scale_y_continuous(expand = expansion(mult = c(0.05, 0.30))) +
   stat_compare_means(
     aes(group = Type),
     method = "wilcox.test",
     paired = TRUE,
     label = "p.signif",
     label.x = 1.45
   ) +
   scale_x_discrete(
     labels = c("A" = "pre", "B" = "post")
   ) +
   theme_bw() +
   theme(legend.position = "none") +
   labs(
     title = "Fraction Cells: pre- vs post-Mepo",
     y = "Fraction of Cells",
     x = ""
   )
 
 library(ggpubr)
 library(dplyr)
 
 # Filter only Macrophages
 cell_fractions_subset <- cell_fractions %>%
   filter(annotation1 %in% c("Macrophages"))
 
 # Check if rows remain
 if (nrow(cell_fractions_subset) == 0) {
   stop("No rows found for Macrophages or Mast cells.")
 }
 
 # Plot
 ggboxplot(cell_fractions_subset, x = "Type", y = "fraction_cells", color = "Type",
   fill = NA,  # Transparent boxes
   add = "jitter",
   facet.by = "annotation1",
   palette = c("grey80", "steelblue"),
   scales = "free_y"
 ) +
   ylab("Cell fraction") +
   stat_compare_means(
     aes(group = Type),
     method = "wilcox.test",
     paired = TRUE,
     label = "p.signif",
     label.x = 1.5
   ) +
   theme(
     panel.border = element_blank(),       # Remove borders
     strip.background = element_blank(),   # Remove facet strip background
     panel.grid = element_blank(),         # Remove grid lines
     panel.background = element_blank(),   # Remove panel background
     axis.line = element_line()            # Keep axes
   )
 
# Initialize a results data frame
results <- data.frame(
  annotation1 = character(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Perform Fisher's exact test for each annotation1
for (label in unique_labels) {
  # Subset data for the current label
  subset_data <- subset(cell_fractions, annotation1 == label)
  
  # Bin fraction_cells into "high" or "low" based on the median
  subset_data$bin <- ifelse(subset_data$fraction_cells > median(subset_data$fraction_cells), "high", "low")
  
  # Create a contingency table between Type and bin
  contingency_table <- table(subset_data$Type, subset_data$bin)
  
  # Perform Fisher's test
  fisher_result <- fisher.test(contingency_table)
  
  # Store results
  results <- rbind(
    results,
    data.frame(annotation1 = label, p_value = fisher_result$p.value, stringsAsFactors = FALSE)
  )
}

# View the results
print(results)

# Optional: Adjust p-values for multiple testing (e.g., FDR correction)
results$adjusted_p_value <- p.adjust(results$p_value, method = "fdr")

# View results with adjusted p-values
print(results)
unique_labels <- unique(cell_fractions$annotation1)

# Initialize a results data frame
results <- data.frame(
  annotation1 = character(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Perform Fisher's exact test for each annotation1
for (label in unique_labels) {
  # Subset data for the current label
  subset_data <- subset(cell_fractions, annotation1 == label)
  
  # Bin fraction_cells into "high" or "low" based on the median
  subset_data$bin <- ifelse(subset_data$fraction_cells > median(subset_data$fraction_cells), "high", "low")
  
  # Create a contingency table between Type and bin
  contingency_table <- table(subset_data$Type, subset_data$bin)
  
  # Perform Fisher's test
  fisher_result <- fisher.test(contingency_table)
  
  # Store results
  results <- rbind(
    results,
    data.frame(annotation1 = label, p_value = fisher_result$p.value, stringsAsFactors = FALSE)
  )
}

# View the results
print(results)

# Optional: Adjust p-values for multiple testing (e.g., FDR correction)
results$adjusted_p_value <- p.adjust(results$p_value, method = "fdr")

# View results with adjusted p-values
print(results)

# Dimplot separating the two conditions
DimPlot(merged_seurat.harmony, reduction = 'umap', group.by = 'annotation1', split.by = "Type", label = TRUE, raster=FALSE)

# Dimplot separating patients and condition
merged_seurat.harmony2 <- merged_seurat.harmony
merged_seurat.harmony2@meta.data$Patient_Type <- paste(merged_seurat.harmony2@meta.data$Patient, merged_seurat.harmony2@meta.data$Type, sep = "_")
DimPlot(merged_seurat.harmony2, reduction = "umap", split.by = "Patient_Type", group.by = 'Patient', raster=FALSE) & NoLegend()

# PSEUDOBULK ANALYSIS
merged_seurat.harmony <- readRDS("~/Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/merged_seurat.harmonyCORRECTED.rds")
merged_seurat.harmony2 <- merged_seurat.harmony
library(DESeq2) # If you have problem with this, ensure you load DESeq2 firs. Save the workspace, close R and reopen the file and load DESeq2
library(ExperimentHub)
library(Seurat)
library(tidyverse)

# pseudo-bulk workflow 
# Acquiring necessary metrics for aggregation across cells in a sample
# 1. counts matrix - sample level
# counts aggregate to sample level

# merged_seurat.harmony2 <- merged_seurat.harmony # Done it above
# View(merged_seurat.harmony2@meta.data)
merged_seurat.harmony2$samples <- paste0(merged_seurat.harmony2$Type, merged_seurat.harmony2$Patient)
merged_seurat.harmony2$cluster_id <- paste0(factor(merged_seurat.harmony2@active.ident))

DefaultAssay(merged_seurat.harmony2)

cts <- AggregateExpression(merged_seurat.harmony2, 
                           group.by = c("annotation1", "samples"),
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
cts.split <- split.data.frame(cts.t, f = factor(splitRows))
# fix colnames and transpose
cts.split.modified <- lapply(cts.split, function(x){
  rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
  t(x)
})

#gsub('.*_(.*)', '\\1', 'B cells_ctrl101')



# Run DE analysis #
# 1. Get counts matrix
counts_cell <- cts.split.modified$`Mast cells`

# 2. generate sample level metadata
colData <- data.frame(samples = colnames(counts_cell))

colData <- colData %>%
  mutate(condition = ifelse(grepl('B', samples), 'PostMEPO', 'PreMEPO')) %>%
  column_to_rownames(var = 'samples')

colData$patient <- substr(rownames(colData), 2, nchar(rownames(colData))) # *

# perform DESeq2 
# Create DESeq2 object   
dds <- DESeqDataSetFromMatrix(countData = counts_cell,
                              colData = colData,
                              design = ~ patient + condition) # try this to paired analysis ~ Patient + condition, needs *

# filter
dds <- dds[rowSums(counts(dds)) > 0,] 
keep <- rowSums(counts(dds) >= 50) >= 12

# or keep <- rowSums(counts(dds)) >=10
#     dds <- dds[keep,]

table(keep)
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "PreMEPO")
# run DESeq2
dds <- DESeq(dds)

summary(results(dds, alpha=0.05))

# Check the coefficients for the comparison
resultsNames(dds)

# Generate results object
res <- results(dds, name = "condition_PostMEPO_vs_PreMEPO")

# save res file
write.csv(res, file = "~/Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/Res-200-5000-paired/basal50h.csv")

# How many adjusted p-values were less than 0.1?
sum(res$padj < 0.05, na.rm=TRUE)
# Histogram pval vs padj
hist(res$pvalue, col = "lavender", border = "white", xlab = "Value", main =  "pval vs padj")
hist(res$padj, col = "steelblue", border = "white", add = TRUE, xlab = "Value", main =  "pval vs pad")
abline( v = 0.05, col = "green4", lwd = 2)
# Plot PCA
rld <- rlog(dds)
plotPCA(rld, intgroup=c("condition")) + geom_text_repel(aes(label=rld@colData$patient))
plotPCA(rld, intgroup=c("patient")) + geom_text_repel(aes(label=rld@colData$patient))

# save res file
# write.csv(res, file = "~/Desktop/scRNA-seq_Brush1/IRIS/res_mast_pre-post.csv")
res$Gene <- rownames(res)
res$neg_log10_padj <- -log10(res$padj)  
res$Significance <- "Not Significant"        

# Mark significant genes based on the conditions
res$Significance[res$padj < 0.05 & res$log2FoldChange > 0] <- "Upregulated"
res$Significance[res$padj < 0.05 & res$log2FoldChange < 0] <- "Downregulated"

colors <- c("Not Significant" = "grey", 
            "Upregulated" = "red", 
            "Downregulated" = "royalblue3")

ggplot(res, aes(x = log2FoldChange, y = neg_log10_padj, color = Significance)) +
  geom_point(size = 2, alpha = 0.6) +  
  scale_color_manual(values = colors) +  
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    title = "DE Ciliated post Mepo"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  
    legend.position = "none"                              
  ) +
  geom_text_repel(
    data = subset(res, padj < 0.05 & abs(log2FoldChange) > 0),  # Subset significant genes
    aes(label = Gene),  
    size = 3,           
    box.padding = 0.3,  
    point.padding = 0.2 
  ) +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed", color = "black") +  # Threshold lines for log2FC
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")  # Threshold line for p-adj

#
#

# Comparison in featureplot
FeaturePlot(merged_seurat.harmony, features = c('IL5RA', 'CSF2RB'), min.cutoff = 'q9', order = TRUE)
FeaturePlot(merged_seurat.harmony, features = c('IL5RA', 'CSF2RB'), split.by = 'Type', min.cutoff = 'q9', order = TRUE)
FeaturePlot(merged_seurat.harmony, features = c('CSF2RB'), split.by = 'Type', min.cutoff = 'q9', order = TRUE)
FeaturePlot(merged_seurat.harmony, features = c("IL4", "IL13", "IL25", "IL33", "TSLP", "IL17RB"), min.cutoff = 'q9', order = TRUE)
FeaturePlot(merged_seurat.harmony, features = c("IL5", "IL6", "CXCL8", "CCL5", "CCL26", "CCL24"), min.cutoff = 'q9', order = TRUE)
FeaturePlot(merged_seurat.harmony, features = c("IFNG", "IL10", "GATA3", "CLCA1", "CCL22", "TGFB1"), min.cutoff = 'q9', order = TRUE)

library("Nebulosa")
library("SCpubr")
plot_density(merged_seurat.harmony, c("IL5RA", 'CSF2RB'), reduction = "umap")
VlnPlot(merged_seurat.harmony, features = c('IL5RA', 'CSF2RB'), group.by = 'annotation1', pt.size = 0)
VlnPlot(merged_seurat.harmony, features = c("VIM", "SMAD2", "IL5RA", "CSF2RB"), group.by = 'annotation1', split.by = 'Type', pt.size = .5, alpha = 0.1)
RidgePlot(merged_seurat.harmony, features = c('IL5RA', 'CSF2RB'), group.by = 'annotation1', ncol = 2)

plot_density(merged_seurat.harmony, c("TRPM5"), reduction = "umap")
VlnPlot(merged_seurat.harmony, features = c("TRPM5"), group.by = 'seurat_clusters', pt.size = 0) # For Tuft cells research

plot_density(merged_seurat.harmony, c("TGFB2"), reduction = "umap") +
facet_grid(.~merged_seurat.harmony@meta.data$Type)

plot_density(merged_seurat.harmony2, c("CLDN5"), reduction = "umap") +
  facet_grid(. ~ merged_seurat.harmony2@meta.data$Type) +
  theme(strip.background = element_rect(fill = "grey", color = "black"),  # Grey box with black border
        strip.text = element_text(color = "black"))  # White text for contrast

#
# Extract % of cells expressing IL5RA per cluster
GetAssayData(merged_seurat.harmony, assay = "RNA", slot = "data")["IL5RA", filtered_cells]

library(Seurat)
library(dplyr)
library(ggplot2)

# Your Seurat object
obj <- merged_seurat.harmony

# Define expression threshold (typical: >0 counts = expressed)
# Define expression threshold
min_reads <- 0
obj$IL5RA_positive <- GetAssayData(merged_seurat.harmony, assay = "RNA", slot = "data")["IL5RA",] > min_reads

# Calculate % IL5RA+ cells per cluster (annotation1)
pct_df <- obj@meta.data %>%
  group_by(annotation1) %>%
  summarise(
    total_cells = n(),
    IL5RA_pos = sum(IL5RA_positive),
    pct_IL5RA = 100 * IL5RA_pos / total_cells
  )

# Order clusters by decreasing percentage
pct_df$annotation1 <- factor(pct_df$annotation1,
                             levels = pct_df$annotation1[order(pct_df$pct_IL5RA, decreasing = TRUE)])

print(pct_df)

# Plot
ggplot(pct_df, aes(x = annotation1, y = pct_IL5RA)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(pct_IL5RA, 1)),   # print values (rounded to 1 decimals)
            vjust = -0.5,                       # move text slightly above the bar
            size = 3) +
  theme_bw() +
  ylab("Percent IL5RA+ cells") +
  xlab("Cluster (annotation1)") +
  ggtitle(paste0("Percent cells expressing IL5RA (min_reads > ", min_reads, ")")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#

# Determine the average expression of a specific gene limited to Basal and Ciliated cells
library(ggplot2)
library(dplyr)

# Subset meta.data to include only Ciliated and Basal cells
filtered_meta <- merged_seurat.harmony@meta.data %>%
  filter(annotation1 %in% c("Ciliated", "Basal"))

# Get the cell names that match the filter
filtered_cells <- rownames(filtered_meta)

# Extract normalized expression for IL5RA only for filtered cells
normalized_counts <- GetAssayData(merged_seurat.harmony, assay = "RNA", slot = "data")["IL5RA", filtered_cells]

# Create a data frame with matched cell metadata
final_df <- data.frame(
  Cell = filtered_cells,
  Expression = as.numeric(normalized_counts),
  Patient = filtered_meta$Patient,   # Add Patient info
  Type = filtered_meta$Type          # Add Type info
)

# Check the first few rows
head(final_df)

# Compute average expression per Patient per Type
avg_expr_df <- final_df %>%
  group_by(Patient, Type) %>%
  summarise(Average_Expression = mean(Expression, na.rm = TRUE), .groups = "drop")

# View the aggregated data
head(avg_expr_df)


p_bar <- ggplot(avg_expr_df, aes(x = Patient, y = Average_Expression, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 1)) +  # Adjust spacing
  theme_minimal() +
  labs(
    title = "Average IL5RA Expression Per Patient and Condition",
    subtitle = "Basal and Ciliated",
    x = "Patient",
    y = "Average Normalized Expression"
  ) +
  scale_fill_manual(
    values = setNames(c("azure3", "coral2"), unique(avg_expr_df$Type))  # Dynamically match values
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 0.5),  # Rotate x-axis labels for clarity
    legend.position = "right"
  )


p_bar

p_box <- ggplot(final_df, aes(x = Patient, y = Expression, fill = Type)) +
  geom_boxplot(
    outlier.shape = NA, 
    position = position_dodge(width = 1),  # Space between boxes
    alpha = 0.2                             # Set box transparency (0 = fully transparent)
  ) +  
  geom_jitter(aes(color = Type), position = position_dodge(width = 1), alpha = 0.5) +  # Align dots with boxes
  stat_summary(
    fun = "mean", 
    geom = "point", 
    aes(color = "black"), 
    position = position_dodge(width = 1),   # Ensure the points align with boxplots
    size = 3,                                  # Size of the mean points
    shape = "-",                                # Shape of the mean point (optional)
    stroke = 2                                  # Border thickness of the mean point
  ) + 
  theme_minimal() +
  labs(
    title = "IL5RA Expression Per Patient and Condition",
    subtitle = "Basal and Ciliated",
    x = "Patient",
    y = "Normalized Expression"
  ) +
  scale_fill_manual(values = setNames(c("azure3", "coral2"), unique(final_df$Type))) +  # Match colors dynamically
  scale_color_manual(values = setNames(c("azure3", "coral2"), unique(final_df$Type))) +  # Match jitter points
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right"
  )

p_box





DotPlot(merged_seurat.harmony, features = c("KRT5", "MKI67", "KRT13", "KRT14", "TP63", "IL33", "KRT8", "SCGB1A1", "MUC5AC", "MUC5B",
                                            "WFDC2", "FOXJ1", "CDHR3", "DNAH11", "CHGA", "ASCL1", "FOOXI1", "CFTR", "CLCNKB",
                                            "POU2F3", "G0S2", "VIM", "CTNNB1", "TOB1", "FOLR1", "DMBT1", "SCGB3A1"), group.by = 'annotation1') + coord_flip()

DotPlot(merged_seurat.harmony, features = c("KRT5", "MKI67", "TP63", "IL33", "KRT8", "SCGB1A1", "MUC5AC", "MUC5B",
                                            "WFDC2", "FOXJ1", "CDHR3", "DNAH11", "FOOXI1", "CFTR", "CLCNKB",
                                            "VIM", "CTNNB1", "TOB1", "SCGB3A1"), group.by = 'annotation1') + coord_flip()




DotPlot(merged_seurat.harmony, c("TP63", "KRT5", "KRT15", 
  "DNAH5", "SPEF2", "PIFO", "FOXJ1", "DNAI1", "CCDC40", "RSPH1", "TEKT1", "TUBB4B", "SNTN",
  "SCGB1A1", "CYP2F1", "CCDC50",
  "IRF8", "DEUP1", "FOXN4", "CDC20B",
  "MUC5AC","MUC5B", "SPDEF",
  "LINC01187", "ATP6V1G3", "FOXI1", "TMPRSS11E", "BSND", "SFTPB", "CFTR",
  "CD68", "CD14", "CSF1R", "MSR1", "CD163", "PPARG", "APOE",
  "TPSAB1", "TPSB2", "HPGDS", "KIT",
  "CD2", "PTPRC"), group.by = "annotation1") + coord_flip()








# Filter only Macrophages and Mast cells
cell_fractions_subset <- cell_fractions %>%
  # filter(annotation1 %in% c("Ciliated", "Basal", "Club","Deuterosomal","Goblet", "Ionocytes","Mast cells", "T-cells", "Macrophages", "Mucociliated", "Dendritic"))
  filter(annotation1 %in% c("Macrophages"))

# Check if rows remain
if (nrow(cell_fractions_subset) == 0) {
  stop("No rows found for Macrophages or Mast cells.")
}

# Plot
ggboxplot(cell_fractions_subset, x = "Type", y = "fraction_cells", color = "Type",
          fill = NA,  # Transparent boxes
          add = "jitter",
          facet.by = "annotation1",
          palette = c("grey80", "dodgerblue2"),
          scales = "free_y"
) +
  ylab("Cell fraction") +
  stat_compare_means(
    aes(group = Type),
    method = "wilcox.test",
    paired = TRUE,
    label = "p.signif", # p.format or p.signif
    label.x = 1.5
  ) +
  theme(
    panel.border = element_blank(),       # Remove borders
    strip.background = element_blank(),   # Remove facet strip background
    panel.grid = element_blank(),         # Remove grid lines
    panel.background = element_blank(),   # Remove panel background
    axis.line = element_line()            # Keep axes
  )


merged_seurat.harmony <- readRDS("~/Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/merged_seurat.harmonyCORRECTED.rds")

library(Seurat)
library(dplyr)
library(ggplot2)

# Your Seurat object
obj <- merged_seurat.harmony

# Define expression threshold (typical: >0 counts = expressed)
# Define expression threshold
min_reads <- 0
obj$IL5RA_positive <- GetAssayData(merged_seurat.harmony, assay = "RNA", slot = "data")["IL5RA",] > 0

# Extract IL5RA raw counts (or use slot = "data" if you want normalized values)
il5ra_counts <- GetAssayData(obj, assay = "RNA", slot = "counts")["IL5RA", ]

# Add counts to metadata for easy grouping
obj$IL5RA_counts <- il5ra_counts

# Load dplyr
library(dplyr)

# Compute averages for TRUE vs FALSE IL5RA_positive cells
avg_df <- obj@meta.data %>%
  group_by(IL5RA_positive) %>%
  summarise(
    avg_IL5RA_counts = mean(IL5RA_counts)
  )

print(avg_df)

# Calculate % IL5RA+ cells per cluster (annotation1)
pct_df <- obj@meta.data %>%
  group_by(annotation1) %>%
  summarise(
    total_cells = n(),
    IL5RA_pos = sum(IL5RA_positive),
    pct_IL5RA = 100 * IL5RA_pos / total_cells
  )

# Order clusters by decreasing percentage
pct_df$annotation1 <- factor(pct_df$annotation1,
                             levels = pct_df$annotation1[order(pct_df$pct_IL5RA, decreasing = TRUE)])

print(pct_df)

# Plot
ggplot(pct_df, aes(x = annotation1, y = pct_IL5RA)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(pct_IL5RA, 1)),   # print values (rounded to 1 decimals)
            vjust = -0.5,                       # move text slightly above the bar
            size = 3) +
  theme_bw() +
  ylab("Percent IL5RA+ cells") +
  xlab("Cluster (annotation1)") +
  ggtitle(paste0("Percent cells expressing IL5RA (min_reads > ", min_reads, ")")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Subset ciliated cells, make a new seurat object with those
# ciliated_seurat <- subset(
#   obj,
#   subset = annotation1 == "Ciliated"
# )

# Subset IL5RA positive
ciliated_IL5RA_positive <- subset(
  obj, # with obj, you have all the cell type > important later
  subset = IL5RA_positive == TRUE
)
# Check average expression of IL5RA
mean(GetAssayData(ciliated_IL5RA_positive, assay = "RNA", slot = "data")["IL5RA", ]) # use slot = data for normalised

# Subset IL5RA negative
ciliated_IL5RA_negative <- subset(
  obj,
  subset = IL5RA_positive == FALSE
)
# Check average expression of IL5RA
mean(GetAssayData(ciliated_IL5RA_negative, assay = "RNA", slot = "data")["IL5RA", ]) # use slot = data for normalised

# DESeq2 IL5RA positive
library(DESeq2) # If you have problem with this, ensure you load DESeq2 firs. Save the workspace, close R and reopen the file and load DESeq2
library(ExperimentHub)
library(Seurat)
library(tidyverse)

# pseudo-bulk workflow 
# Acquiring necessary metrics for aggregation across cells in a sample
# 1. counts matrix - sample level
# counts aggregate to sample level

# merged_seurat.harmony2 <- merged_seurat.harmony # Done it above
# View(merged_seurat.harmony2@meta.data)
ciliated_IL5RA_positive$samples <- paste0(ciliated_IL5RA_positive$Type, ciliated_IL5RA_positive$Patient)
ciliated_IL5RA_positive$cluster_id <- paste0(factor(ciliated_IL5RA_positive@active.ident))

DefaultAssay(ciliated_IL5RA_positive)

cts <- AggregateExpression(ciliated_IL5RA_positive, 
                           group.by = c("annotation1", "samples"),
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
cts.split <- split.data.frame(cts.t, f = factor(splitRows))
# fix colnames and transpose
cts.split.modified <- lapply(cts.split, function(x){
  rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
  t(x)
})

#gsub('.*_(.*)', '\\1', 'B cells_ctrl101')



# Run DE analysis #
# 1. Get counts matrix
counts_cell <- cts.split.modified$`Ciliated`

# 2. generate sample level metadata
colData <- data.frame(samples = colnames(counts_cell))

colData <- colData %>%
  mutate(condition = ifelse(grepl('B', samples), 'PostMEPO', 'PreMEPO')) %>%
  column_to_rownames(var = 'samples')

colData$patient <- substr(rownames(colData), 2, nchar(rownames(colData))) # *

# perform DESeq2 
# Create DESeq2 object   
dds <- DESeqDataSetFromMatrix(countData = counts_cell,
                              colData = colData,
                              design = ~ patient + condition) # try this to paired analysis ~ Patient + condition, needs *

# filter
dds <- dds[rowSums(counts(dds)) > 0,] 
keep <- rowSums(counts(dds) >= 10) >= 12

# or keep <- rowSums(counts(dds)) >=10
#     dds <- dds[keep,]

table(keep)
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "PreMEPO")
# run DESeq2
dds <- DESeq(dds)

summary(results(dds, alpha=0.05))

# Check the coefficients for the comparison
resultsNames(dds)

# Generate results object
res <- results(dds, name = "condition_PostMEPO_vs_PreMEPO")

res_positive_10 <- res

write.csv(res_positive_10, "~/Desktop/res_ciliated_positive_10.csv")







# For the IL5RA negative
ciliated_IL5RA_negative$samples <- paste0(ciliated_IL5RA_negative$Type, ciliated_IL5RA_negative$Patient)
ciliated_IL5RA_negative$cluster_id <- paste0(factor(ciliated_IL5RA_negative@active.ident))

DefaultAssay(ciliated_IL5RA_negative)

cts <- AggregateExpression(ciliated_IL5RA_negative, 
                           group.by = c("annotation1", "samples"),
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
cts.split <- split.data.frame(cts.t, f = factor(splitRows))
# fix colnames and transpose
cts.split.modified <- lapply(cts.split, function(x){
  rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
  t(x)
})

#gsub('.*_(.*)', '\\1', 'B cells_ctrl101')



# Run DE analysis #
# 1. Get counts matrix
counts_cell <- cts.split.modified$`Ciliated`

# 2. generate sample level metadata
colData <- data.frame(samples = colnames(counts_cell))

colData <- colData %>%
  mutate(condition = ifelse(grepl('B', samples), 'PostMEPO', 'PreMEPO')) %>%
  column_to_rownames(var = 'samples')

colData$patient <- substr(rownames(colData), 2, nchar(rownames(colData))) # *

# perform DESeq2 
# Create DESeq2 object   
dds <- DESeqDataSetFromMatrix(countData = counts_cell,
                              colData = colData,
                              design = ~ patient + condition) # try this to paired analysis ~ Patient + condition, needs *

# filter
dds <- dds[rowSums(counts(dds)) > 0,] 
keep <- rowSums(counts(dds) >= 10) >= 12

# or keep <- rowSums(counts(dds)) >=10
#     dds <- dds[keep,]

table(keep)
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "PreMEPO")
# run DESeq2
dds <- DESeq(dds)

summary(results(dds, alpha=0.05))

# Check the coefficients for the comparison
resultsNames(dds)

# Generate results object
res <- results(dds, name = "condition_PostMEPO_vs_PreMEPO")

res_negative_10 <- res

write.csv(res_negative_10, "~/Desktop/res_ciliated_negative_10.csv")









res_positive_10 <- read.csv("~/Desktop/res_ciliated_positive_10.csv", header = T)
res_negative_10 <- read.csv("~/Desktop/res_ciliated_negative_10.csv", header = T)

names(res_positive_10)[1] <- "Gene"
names(res_negative_10)[1] <- "Gene"
# total RNA volcanoplot #
library(ggplot2)
library(ggrepel)



res_positive_10$neg_log10_padj <- -log10(res_positive_10$padj)  
res_positive_10$Significance <- "Not Significant"        

# Mark significant genes based on the conditions
res_positive_10$Significance[res_positive_10$padj < 0.05 & res_positive_10$log2FoldChange > 0.1] <- "Upregulated"
res_positive_10$Significance[res_positive_10$padj < 0.05 & res_positive_10$log2FoldChange < -0.1] <- "Downregulated"

colors <- c("Not Significant" = "grey", 
            "Upregulated" = "red", 
            "Downregulated" = "royalblue3")

ggplot(res_positive_10, aes(x = log2FoldChange, y = neg_log10_padj, color = Significance)) +
  geom_point(size = 2, alpha = 0.6) +  
  scale_color_manual(values = colors) +  
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    title = "6 months IL5RA-positive-ciliated"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  
    legend.position = "none",
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  ) +
  geom_text_repel(
    data = subset(res_positive_10, padj < 0.05 & abs(log2FoldChange) > 1),  # Subset significant genes
    aes(label = Gene),  
    size = 3,           
    box.padding = 0.3,  
    point.padding = 0.2 
  ) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "black") +  # Threshold lines for log2FC
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") + # Threshold line for p-adj
  xlim(-3,3)+
  ylim(0, 5)
#
#


# MA plot res_positive_10 RNA #
# Create a column for significance
res_positive_10$Significance <- ifelse(
  res_positive_10$padj < 0.05 & res_positive_10$log2FoldChange > 0.1, "Upregulated",
  ifelse(res_positive_10$padj < 0.05 & res_positive_10$log2FoldChange < -0.1, "Downregulated", "Not significant")
)

ma_colors <- c("Upregulated" = "red", "Downregulated" = "royalblue3", "Not significant" = "grey")

ggplot(res_positive_10, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(data = subset(res_positive_10, Significance == "Not significant"),
             aes(shape = Significance), color = "grey", size = 2, alpha = 0.6) +  
  geom_point(data = subset(res_positive_10, Significance != "Not significant"),
             aes(color = Significance), size = 2, alpha = 0.6) +  
  scale_x_log10() +  # Log scale for baseMean
  scale_color_manual(values = ma_colors) +  
  scale_shape_manual(values = c("Not significant" = 4)) +  
  labs(
    x = "Mean Expression (log10 scale)",
    y = "Log2 Fold Change",
    title = "6 months IL5RA-positive-ciliated"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  
    legend.title = element_blank(),
    legend.position = "none"  
  ) +
  geom_hline(yintercept = c(-0.1, 0, 0.1), linetype = "dashed", color = "black")  

#
#


