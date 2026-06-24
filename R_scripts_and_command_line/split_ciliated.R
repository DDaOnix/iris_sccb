
library(Seurat)
library(dplyr)
library(ggplot2)

merged_seurat.harmony <- readRDS("~/Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/merged_seurat.harmonyCORRECTED.rds")
# Your Seurat object
obj <- merged_seurat.harmony

# Define expression threshold (typical: >0 counts = expressed)
# Define expression threshold
min_reads <- 0
obj$IL5RA_positive <- GetAssayData(obj, assay = "RNA", layer = "data")["IL5RA",] > 0

obj$IL5RA_positive <- ifelse(obj$IL5RA_positive, "T", "F")


# Extract IL5RA raw counts (or use slot = "layer" if you want normalized values)
il5ra_counts <- GetAssayData(obj, assay = "RNA", layer = "counts")["IL5RA", ]


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
# Check number of Ciliated cells
sum(ciliated_IL5RA_positive$annotation1 == "Ciliated")

# Subset IL5RA negative
ciliated_IL5RA_negative <- subset(
  obj,
  subset = IL5RA_positive == FALSE
)
# Check average expression of IL5RA
mean(GetAssayData(ciliated_IL5RA_negative, assay = "RNA", slot = "data")["IL5RA", ]) # use slot = data for normalised
# Check number of Ciliated cells
sum(ciliated_IL5RA_negative$annotation1 == "Ciliated")
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







# Perform DESeq2 on ciliated cells positive vs negative for the expression on IL5RA
library(Seurat)
library(dplyr)
library(ggplot2)

merged_seurat.harmony <- readRDS("~/Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/merged_seurat.harmonyCORRECTED.rds")
# Your Seurat object
obj <- merged_seurat.harmony

# Define expression threshold (typical: >0 counts = expressed)
# Define expression threshold
min_reads <- 0
obj$IL5RA_positive <- GetAssayData(merged_seurat.harmony, assay = "RNA", layer = "data")["IL5RA",] > 0

obj$IL5RA_positive <- ifelse(obj$IL5RA_positive, "T", "F")

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
obj$samples <- paste0(obj$IL5RA_positive, obj$Patient)
obj$cluster_id <- paste0(factor(obj@active.ident))

DefaultAssay(obj)

cts <- AggregateExpression(obj, 
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
  mutate(condition = ifelse(grepl('T', samples), 'Positive', 'Negative')) %>%
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
dds$condition <- relevel(dds$condition, ref = "Negative")
# run DESeq2
dds <- DESeq(dds)

summary(results(dds, alpha=0.05))

# Check the coefficients for the comparison
resultsNames(dds)

# Generate results object
res <- results(dds, name = "condition_Positive_vs_Negative")

res_positive_10 <- res


write.csv(res_positive_10, "~/Desktop/res_ciliated_positive_vs_negative.csv")


# volcanoplot #
library(ggplot2)
library(ggrepel)

res_positive_10 <- read.csv("~/Desktop/res_ciliated_positive_vs_negative.csv", header = T)

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
    title = "Positive vs Negative ciliated for IL5RA"
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
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") #+ # Threshold line for p-adj
# xlim(-3,3)+
# ylim(0, 5)

# MA Plot
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
    title = "Positive vs Negative ciliatd for IL5RA"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  
    legend.title = element_blank(),
    legend.position = "none"  
  ) +
  geom_hline(yintercept = c(-0.1, 0, 0.1), linetype = "dashed", color = "black")  





















######################################################################
# Split A and B from Type and then check the DEGs between +ve and -ve
library(Seurat)
library(dplyr)
library(ggplot2)

merged_seurat.harmony <- readRDS("~/Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/merged_seurat.harmonyCORRECTED.rds")
# Your Seurat object
obj <- merged_seurat.harmony

# Define expression threshold (typical: >0 counts = expressed)
# Define expression threshold
min_reads <- 0
obj$IL5RA_positive <- GetAssayData(merged_seurat.harmony, assay = "RNA", layer = "data")["IL5RA",] > 0

obj$IL5RA_positive <- ifelse(obj$IL5RA_positive, "T", "F")


# Extract IL5RA raw counts (or use slot = "layer" if you want normalized values)
il5ra_counts <- GetAssayData(obj, assay = "RNA", layer = "counts")["IL5RA", ]


# Bdd counts to metadata for easy grouping
obj$IL5RA_counts <- il5ra_counts
# Subset Seurat object by metadata column "Type"

obj_A <- subset(obj, subset = Type == "A")
obj_B <- subset(obj, subset = Type == "B")

# Checks
table(obj$Type)
table(obj_A$Type)
table(obj_B$Type)

# DEGs positive ciliated vs negative in A
obj_A$samples <- paste0(obj_A$IL5RA_positive, obj_A$Patient)
obj_A$cluster_id <- paste0(factor(obj_A@active.ident))

DefaultAssay(obj_A)

cts <- AggregateExpression(obj_A, 
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

#gsub('.*_(.*)', '\\1', 'A cells_ctrl101')

# Run DE analysis #
# 1. Get counts matrix
counts_cell <- cts.split.modified$`Ciliated`

# 2. generate sample level metadata
colData <- data.frame(samples = colnames(counts_cell))

colData <- colData %>%
  mutate(condition = ifelse(grepl('T', samples), 'Positive', 'Negative')) %>%
  column_to_rownames(var = 'samples')

colData$patient <- substr(rownames(colData), 2, nchar(rownames(colData))) # *

# perform DESeq2 
# Create DESeq2 obj_Aect   
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
dds$condition <- relevel(dds$condition, ref = "Negative")
# run DESeq2
dds <- DESeq(dds)

summary(results(dds, alpha=0.05))

# Check the coefficients for the comparison
resultsNames(dds)

# Generate results obj_Aect
res <- results(dds, name = "condition_Positive_vs_Negative")
# save res
write.csv(res, "~/Desktop/res_ciliated_positive_vs_negative_A.csv")


# DEGs positive ciliated vs negative in B
obj_B$samples <- paste0(obj_B$IL5RA_positive, obj_B$Patient)
obj_B$cluster_id <- paste0(factor(obj_B@active.ident))

DefaultAssay(obj_B)

cts <- AggregateExpression(obj_B, 
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
  mutate(condition = ifelse(grepl('T', samples), 'Positive', 'Negative')) %>%
  column_to_rownames(var = 'samples')

colData$patient <- substr(rownames(colData), 2, nchar(rownames(colData))) # *

# perform DESeq2 
# Create DESeq2 obj_Bect   
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
dds$condition <- relevel(dds$condition, ref = "Negative")
# run DESeq2
dds <- DESeq(dds)

summary(results(dds, alpha=0.05))

# Check the coefficients for the comparison
resultsNames(dds)

# Generate results obj_B
res <- results(dds, name = "condition_Positive_vs_Negative")
# save res
write.csv(res, "~/Desktop/res_ciliated_positive_vs_negative_B.csv")

# Load res_A and res_B
res_A <- read.csv("~/Desktop/res_ciliated_positive_vs_negative_A.csv", header = T)
res_B <- read.csv("~/Desktop/res_ciliated_positive_vs_negative_B.csv", header = T)

names(res_A)[1] <- "Gene"
names(res_B)[1] <- "Gene"


library(ggplot2)
library(ggrepel)
# Volcanoplot A
res_A$neg_log10_padj <- -log10(res_A$padj)  
res_A$Significance <- "Not Significant"        

# Mark significant genes based on the conditions
res_A$Significance[res_A$padj < 0.05 & res_A$log2FoldChange > 0.1] <- "Upregulated"
res_A$Significance[res_A$padj < 0.05 & res_A$log2FoldChange < -0.1] <- "Downregulated"

colors <- c("Not Significant" = "grey", 
            "Upregulated" = "red", 
            "Downregulated" = "royalblue3")

ggplot(res_A, aes(x = log2FoldChange, y = neg_log10_padj, color = Significance)) +
  geom_point(size = 2, alpha = 0.6) +  
  scale_color_manual(values = colors) +  
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    title = "Positive vs Negative ciliated for IL5RA in A"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  
    legend.position = "none",
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  ) +
  geom_text_repel(
    data = subset(res_A, padj < 0.05 & abs(log2FoldChange) > 1),  # Subset significant genes
    aes(label = Gene),  
    size = 3,           
    box.padding = 0.3,  
    point.padding = 0.2 
  ) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "black") +  # Threshold lines for log2FC
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") 



# Volcanoplot B
res_B$neg_log10_padj <- -log10(res_B$padj)  
res_B$Significance <- "Not Significant"        

# Mark significant genes based on the conditions
res_B$Significance[res_B$padj < 0.05 & res_B$log2FoldChange > 0.1] <- "Upregulated"
res_B$Significance[res_B$padj < 0.05 & res_B$log2FoldChange < -0.1] <- "Downregulated"

colors <- c("Not Significant" = "grey", 
            "Upregulated" = "red", 
            "Downregulated" = "royalblue3")

ggplot(res_B, aes(x = log2FoldChange, y = neg_log10_padj, color = Significance)) +
  geom_point(size = 2, alpha = 0.6) +  
  scale_color_manual(values = colors) +  
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 Bdjusted P-value",
    title = "Positive vs Negative ciliated for IL5RA in B"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  
    legend.position = "none",
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  ) +
  geom_text_repel(
    data = subset(res_B, padj < 0.05 & abs(log2FoldChange) > 1),  # Subset significant genes
    aes(label = Gene),  
    size = 3,           
    box.padding = 0.3,  
    point.padding = 0.2 
  ) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "black") +  # Threshold lines for log2FC
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") 




















# DESeq2 with multiple parameters to compare: B:A, positive:negative

library(Seurat)
library(dplyr)
library(ggplot2)

obj <- readRDS("~/Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/merged_seurat.harmonyCORRECTED.rds")

# Define expression threshold (typical: >0 counts = expressed)

min_reads <- 0
obj$IL5RA_positive <- GetAssayData(obj, assay = "RNA", layer = "data")["IL5RA",] > 0

obj$IL5RA_positive <- ifelse(obj$IL5RA_positive, "T", "F")


# Extract IL5RA raw counts (or use slot = "layer" if you want normalized values)
il5ra_counts <- GetAssayData(obj, assay = "RNA", layer = "counts")["IL5RA", ]


# Add counts to metadata for easy grouping
obj$IL5RA_counts <- il5ra_counts

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
obj$samples <- paste0(obj$IL5RA_positive, obj$Type, obj$Patient)
obj$cluster_id <- paste0(factor(obj@active.ident))

DefaultAssay(obj)

cts <- AggregateExpression(obj, 
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



# Run DE analysis #
# 1. Get counts matrix
counts_cell <- cts.split.modified$`Ciliated`

# 2. generate sample level metadata
colData <- data.frame(samples = colnames(counts_cell))

colData <- colData %>%
  mutate(genotype = ifelse(grepl('T', samples), 'Positive', 'Negative')) %>%
  mutate(condition = ifelse(grepl('B', samples), 'post_Mepo', 'pre_Mepo')) %>%
  column_to_rownames(var = 'samples')
colData$patient <- sub(".*P", "P", rownames(colData)) # Take from the rownames everyting after P, P included and create a new column called patient


# perform DESeq2 
# Create DESeq2 object   
dds <- DESeqDataSetFromMatrix(countData = counts_cell,
                              colData = colData,
                              design = ~ patient + condition + genotype +condition:genotype) # try this to paired analysis ~ Patient + genotype, needs *

dds$condition <- relevel(dds$condition, ref = "pre_Mepo")
dds$genotype <- relevel(dds$genotype, ref = "Negative")
# filter
dds <- dds[rowSums(counts(dds)) > 0,] 
keep <- rowSums(counts(dds) >= 10) >= 12

# or keep <- rowSums(counts(dds)) >=10
#     dds <- dds[keep,]

table(keep)
dds <- dds[keep,]

# run DESeq2
dds <- DESeq(dds)

summary(results(dds, alpha=0.05))

# Check the coefficients for the comparison
resultsNames(dds)

results(dds, contrast = c("condition", "post_Mepo", "pre_Mepo"))
results(dds, contrast = c("genotype", "Positive", "Negative"))
results(dds, contrast = list(c("condition_post_Mepo_vs_pre_Mepo", "conditionpost_Mepo.genotypePositive")))
