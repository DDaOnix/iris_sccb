# Try scplotter scplotter website
# scplotter is an R package that is built upon plotthis. It provides a set of functions to visualize single-cell sequencing data in an easy and efficient way.
# https://github.com/pwwang/scplotter/?tab=readme-ov-file

# Installation 
remotes::install_github("pwwang/scplotter")
# or
devtools::install_github("pwwang/scplotter")
# or using conda
 conda install pwwang::r-scplotter


library(scplotter)



rows_to_label <- c("HLA-DRA", "HLA-DPB1", "HLA-DQB1", "HLA-DPA1",  "HLA-DRB1",  "HLA-DMB")  # Replace with your actual row names

data_1 = b.mepo.response[rownames(b.mepo.response) %in% rows_to_label, ]

ggplot(b.mepo.response, aes(x = avg_log2FC, y = neg_log10_pvalue)) +
  geom_point(aes(color = b.mepo.response[rownames(b.mepo.response) %in% rows_to_label, ]), size = 2, alpha = 0.6) +  # Highlight significant points
  scale_color_manual(values = c( "red")) +  # Use red for significant points
  labs(
    x = "Log2 Fold Change",
    y = "-Log10(P-value)",
    title = "Basal resting Pre vs Post Mepo"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +  # Add fold change threshold lines
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # Add p-value threshold line
  geom_text_repel(
    data_1 = b.mepo.response[rownames(b.mepo.response) %in% rows_to_label, ],  # Subset rows by rownames
    aes(label = rownames(b.mepo.response)[rownames(b.mepo.response) %in% rows_to_label]),
    segment.color = 'black',
    size = 3,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = Inf
  )









library(ggplot2)
library(ggrepel)

# Define row names to be labeled
rows_to_label <- c("HLA-DRA", "HLA-DPB1", "HLA-DQB1", "HLA-DPA1",  "HLA-DRB1",  "HLA-DMB", "HLA-DRB5")  # Replace with your actual row names




# Add a new column 'is_labeled' in b.mepo.response to identify points to color
b.mepo.response2 <- b.mepo.response
b.mepo.response2$is_labeled <- ifelse(rownames(b.mepo.response2) %in% rows_to_label, "highlight", "normal")

ggplot(b.mepo.response2, aes(x = avg_log2FC, y = neg_log10_pvalue)) +
  geom_point(aes(color = is_labeled), size = 2, alpha = 0.6) +  # Use 'is_labeled' to color points
  scale_color_manual(values = c("normal" = "grey90", "highlight" = "red")) +  # Color only highlighted points red
  labs(
    x = "Log2 Fold Change",
    y = "-Log10(P-value)",
    title = "Basal resting Pre vs Post Mepo"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +  # Add fold change threshold lines
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # Add p-value threshold line
  geom_text_repel(
    data = b.mepo.response[rownames(b.mepo.response) %in% rows_to_label, ],  # Subset rows to label
    aes(label = rownames(b.mepo.response)[rownames(b.mepo.response) %in% rows_to_label]),
    segment.color = 'transparent',
    size = 3,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = Inf
  )











# pseudobulk the counts based on donor-condition-celltype
pseudo_iris <- AggregateExpression(merged_seurat.harmony2, assays = "RNA", return.seurat = T, group.by = c("Type", "Patient", "clustering3"))

# each 'cell' is a donor-condition-celltype pseudobulk profile
tail(Cells(pseudo_iris))

pseudo_iris$celltype.stim <- paste(pseudo_iris$clustering3, pseudo_iris$Type, sep = "_")

Idents(pseudo_iris) <- "celltype.stim"

bulk.basal.de <- FindMarkers(object = pseudo_iris, 
                            ident.1 = "Mast cells_B", 
                            ident.2 = "Mast cells_A",
                            test.use = "MAST")
head(bulk.basal.de, n = 15)

# Load ggplot2
library(ggplot2)
library(ggrepel)
# Transform p-value column to -log10 scale for the y-axis
bulk.basal.de$neg_log10_pvalue <- -log10(bulk.basal.de$p_val_adj)

# Create a new column to identify significant genes based on p-value and fold change thresholds
bulk.basal.de$significant <- bulk.basal.de$p_val_adj < 0.05 & abs(bulk.basal.de$avg_log2FC) > 1


# Save the table
#write.csv(bulk.basal.de, "/Volumes/prj_id_iris/IRIS/Healthy_ref_GSE143868/de_basal_resting_pre_vs_post.csv")

data = subset(bulk.basal.de, significant)
# rows_to_label <- c("HLA-DRA", "HLA-DPB1", "HLA-DQB1", "HLA-DPA1",  "HLA-DRB1",  "HLA-DMB")  # Replace with your actual row names
# bulk.basal.de2$is_labeled <- ifelse(rownames(bulk.basal.de2) %in% rows_to_label, "highlight", "normal") # If you want to colour only certain genes in rows_to_label
# Create the volcano plot
ggplot(bulk.basal.de, aes(x = avg_log2FC, y = neg_log10_pvalue)) +
  geom_point(aes(color = (p_val_adj < 0.05 & abs(avg_log2FC) > 1)), size = 2, alpha = 0.6) +  # Highlight significant points | colour = is_labelled to colour only selected genes
  scale_color_manual(values = c("grey", "red")) +  # Use red for significant points | values = c("normal" = "grey90", "highlight" = "red")
  labs(
    x = "Log2 Fold Change",
    y = "-Log10(P-value)",
    title = "Mast Pre vs Post Mepo"
  ) +
  #xlim(-4, 4) +  # Set x-axis limits
  # ylim(0, 3000) +  # Set y-axis limits
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"  
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +  # Add fold change threshold lines
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") # +  # Add p-value threshold line
  geom_text_repel(segment.color = 'transparent',

                  data = subset(bulk.basal.de, significant),
                  aes(label = rownames(data)), # Label only significant genes

                  # data = bulk.basal.de[rownames(bulk.basal.de) %in% rows_to_label, ],  # Subset rows by rownames
                  # aes(label = rownames(bulk.basal.de)[rownames(bulk.basal.de) %in% rows_to_label]),

                  size = 3,
                  box.padding = 0.3,
                  point.padding = 0.2,
                  max.overlaps = Inf
  )
##


install.packages("scCustomize")
install.packages("BiocManager")
BiocManager::install("Nebulosa")
library(scCustomize)

scCustomize::Plot_Density_Custom(seurat_object = merged_seurat.harmony, features = "CD3D",
                                 viridis_palette= "viridis")




library(Seurat)

# Specify the replicate you want to visualize
selected_replicate <- "P18"  # Replace with the specific replicate name

# Subset the Seurat object to include only cells from the selected replicate
subset_seurat <- subset(merged_seurat.harmony, subset = Patient == selected_replicate)

# Generate the DimPlot, splitting by condition
DimPlot(subset_seurat, reduction = "umap", split.by = "Type") + 
  ggtitle(paste("DimPlot for", selected_replicate))  # Optional title






ggplot(cell_fractions, aes(x = Type, y = fraction_cells, fill = annotation1)) + # modify accordingly to the df used cell_counts or cell_fractions
  geom_boxplot() +
  facet_wrap(~ annotation1, scales = "free_y", nrow = 1) +  # Create a separate box plot for each cell type
  labs(
    x = "Condition",
    y = "Fraction cells", # Modify accordingly Cell number of Cell fraction
    title = "Fraction of Cells per Cell Type Across Conditions" # Modify accordingly Cell number of Cell fraction
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  ) +
  #scale_fill_brewer(palette = "Greys") +  # Optional: Use a color palette for conditions
  stat_compare_means(method = "wilcox.test",    # Specify the statistical test (e.g., "t.test", "wilcox.test")
                     label = "p.signif",     # Use "p.signif" for significance stars or "p.format" for p-values
                     comparisons = list(c("A", "B"))  # Pairwise comparisons
  )

unique_labels <- unique(cell_fractions$annotation1)



















# PSEUDOBULK ANALYSIS
# pseudo-bulk workflow -----------------
# Acquiring necessary metrics for aggregation across cells in a sample
# 1. counts matrix - sample level
# counts aggregate to sample level
merged_seurat.harmony2 <- merged_seurat.harmony

View(merged_seurat.harmony2@meta.data)

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
counts_cell <- cts.split.modified$`Basal`

# 2. generate sample level metadata
colData <- data.frame(samples = colnames(counts_cell))

colData <- colData %>%
  mutate(condition = ifelse(grepl('B', samples), 'PostMEPO', 'PreMEPO')) %>%
  column_to_rownames(var = 'samples')
#

  colData$patient_id <- substring(rownames(colData), 2)


# get more information from metadata


# perform DESeq2 --------
# Create DESeq2 object   
dds <- DESeqDataSetFromMatrix(countData = counts_cell,
                              colData = colData,
                              design = ~ patient_id + condition) # try this to paired analysis ~ Patient + condition


# dds <- DESeqDataSetFromMatrix(cluster_counts,colData = cluster_metadata,design = ~ patient_id + treatment) # Suggested by Rocio to force the Paired analysis


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
sum(res$padj < 0.05, na.rm=TRUE)

# save res file
# write.csv(res, file = "~/Desktop/scRNA-seq_Brush1/IRIS/res_mast_pre-post.csv")

# generate volcanoplot
library(EnhancedVolcano)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                #xlim = c(-5, 5),
                title = "DE Basal pre/post MEPO corrected",
                pCutoff = 0.05,
                FCcutoff = 1,
                col = c("grey", "grey", 'grey', "red3"),
                #selectLab = c("GSTA1", "GSTA2","HLA-DQA1", "HLA-DRB1", "NQO1", "B2M", "PRDX1"),
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                drawConnectors = FALSE,
                boxedLabels = FALSE,
                legendPosition = "top")









pseudo_iris$celltype.stim <- paste(pseudo_iris$annotation1, pseudo_iris$Type, sep = "_")

Idents(pseudo_iris) <- "celltype.stim"

bulk.cell_type.de <- FindMarkers(object = pseudo_iris, 
                                 ident.1 = "Basal_B", 
                                 ident.2 = "Basal_A",
                                 test.use = "DESeq2")
head(bulk.cell_type.de, n = 15)

# Load ggplot2
library(ggplot2)
library(ggrepel)
# Transform p-value column to -log10 scale for the y-axis
bulk.cell_type.de$neg_log10_pvalue <- -log10(bulk.cell_type.de$p_val)

# Create a new column to identify significant genes based on p-value and fold change thresholds
bulk.cell_type.de$significant <- bulk.cell_type.de$p_val < 0.05 & abs(bulk.cell_type.de$avg_log2FC) > 1


# Save the table
#write.csv(bulk.cell_type.de, "/Volumes/prj_id_iris/IRIS/Healthy_ref_GSE143868/de_basal_resting_pre_vs_post.csv")

data = subset(bulk.cell_type.de, significant)
# rows_to_label <- c("HLA-DRA", "HLA-DPB1", "HLA-DQB1", "HLA-DPA1",  "HLA-DRB1",  "HLA-DMB")  # Replace with your actual row names
# bulk.cell_type.de2$is_labeled <- ifelse(rownames(bulk.cell_type.de2) %in% rows_to_label, "highlight", "normal") # If you want to colour only certain genes in rows_to_label
# Create the volcano plot
ggplot(bulk.cell_type.de, aes(x = avg_log2FC, y = neg_log10_pvalue)) +
  geom_point(aes(color = (p_val_adj < 0.05 & abs(avg_log2FC) > 1)), size = 2, alpha = 0.6) +  # Highlight significant points | colour = is_labelled to colour only selected genes
  scale_color_manual(values = c("grey", "red")) +  # Use red for significant points | values = c("normal" = "grey90", "highlight" = "red")
  labs(
    x = "Log2 Fold Change",
    y = "-Log10(P-adj)",
    title = "TEST DE Post Mepo"
  ) +
  #xlim(-4, 4) +  # Set x-axis limits
  # ylim(0, 3000) +  # Set y-axis limits
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"  
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +  # Add fold change threshold lines
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") # +  # Add p-value threshold line
geom_text_repel(segment.color = 'transparent',
                
                data = subset(bulk.cell_type.de, significant),
                aes(label = rownames(data)), # Label only significant genes
                
                data = bulk.cell_type.de[rownames(bulk.cell_type.de) %in% rows_to_label, ],  # Subset rows by rownames
                aes(label = rownames(bulk.cell_type.de)[rownames(bulk.cell_type.de) %in% rows_to_label]),
                
                size = 3,
                box.padding = 0.3,
                point.padding = 0.2,
                max.overlaps = Inf
)
