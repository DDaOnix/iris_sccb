# Compare DE genes IL-5 stimulated and Ciliated Mepo treated 

#Load tables
cilDE <- read.csv("~/Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/Res-200-5000-paired/ciliated50h.csv", header = T)
IL5t <- read.csv("~/Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/IL_5_resp_significant.csv", header = T) # From https://doi.org/10.1111/all.14297

# Modify the tables to get just two columns: Gene and log2FoldChange_IL5 or log2FoldChange_MepoCil
names(cilDE)[1] <- "Gene" # Name the first column
names(cilDE)[3] <- "log2FoldChange_MepoCil" # Rename the column Foldchange
cilDE <- cilDE[,-c(2,4,5,6,7,8)] # Remove other columns
names(IL5t)[2] <- "log2FoldChange_IL5"
IL5t <- IL5t[, -c(1,3,4,5)]

# This is a list of genes found to be DE in the table IL5 stimulation and in the ciliated cells from IRIS study
common_genes <- data.frame(
  Gene = genes <- c(
    "ABI1", "ACVR2A", "ADGRA3", "AKR1A1", "ALCAM", "ANKRD17", "ANOS1", "AP1S1",
    "ARHGAP21", "ARHGAP42", "ARHGEF12", "ATP6V0D1", "BCL2L11", "BMP2K", "CCDC28A",
    "CCDC89", "CCT3", "CCT5", "CDKAL1", "CDYL", "CENPP", "CFL1", "CHD7", "CLDN16",
    "CREB3", "CXXC5", "DDX10", "DENND1B", "DIAPH2", "DMTF1", "DUS1L", "DYRK1A",
    "ERCC8", "FAR2", "FTH1", "FUCA2", "GADD45A", "GANC", "GAPDH", "GDA", "GLB1L",
    "GPATCH8", "HMGCS1", "HNF4G", "HTATIP2", "KRT8", "LDLR", "LRCH3", "LRP10",
    "MAGI3", "MAPK8", "MBD5", "MBTD1", "MID1IP1", "MMS19", "MRPS6", "MYO5A",
    "NCOA1", "NR3C1", "NUP160", "NXN", "ORC5", "PCBD2", "PDCL3", "PDLIM1", "PGAM1",
    "PHLPP1", "PLAG1", "PPARGC1A", "RAB11FIP2", "ROBO1", "RPRD1A", "SCD5", "SH3RF1",
    "SIK3", "SMG7", "SNRPA1", "ST8SIA4", "STAG1", "STAU1", "STK3", "SULT1A1",
    "TEAD1", "THADA", "THRB", "TPM2", "TPP2", "TRIM2", "TYW1B", "UBB", "USP15",
    "USP24", "USP31", "USP32", "VAMP3", "WDR7", "ZCCHC7", "ZDHHC21"
  )
  
)

# Merge the tables to include only the genes in common
merged_table <- merge(common_genes, IL5t, by = "Gene")  # Add FoldChange1
merged_table <- merge(merged_table, cilDE, by = "Gene")  # Add FoldChange2


library(pheatmap)
# Prepare the data for the heat map
# Assuming `merged_table` contains columns: "Gene", "log2FoldChange_IL5", "log2FoldChange_MepoCil"
# Set row names as genes and select the two fold change columns
heatmap_data <- merged_table[, c("log2FoldChange_IL5", "log2FoldChange_MepoCil")]
rownames(heatmap_data) <- merged_table$Gene

# Generate the heat map
# Determine the limits for the color scale
heatmap_range <- max(abs(heatmap_data), na.rm = TRUE)

# pheatmap(
#   mat = heatmap_data,
#   scale = "none",
#   cluster_rows = TRUE,
#   cluster_cols = TRUE,
#   color = colorRampPalette(c("blue", "white", "red"))(100),
#   breaks = seq(-heatmap_range, heatmap_range, length.out = 101),
#   main = "DE IL-5 stim vs. Mepo (Ciliated cells)",
#   angle_col = 0,
#   border_color = NA,
#   fontsize_row = 6
# )
# Cap the heatmap data at -1 and 1
heatmap_data_capped <- heatmap_data
heatmap_data_capped[heatmap_data_capped < -1] <- -1
heatmap_data_capped[heatmap_data_capped > 1] <- 1

# Define color and breaks
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- seq(-1, 1, length.out = 101)

# Draw the heatmap
pheatmap(
  mat = heatmap_data_capped,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = color_palette,
  breaks = breaks,
  main = "DE IL-5 stim vs. Mepo (Ciliated cells)",
  angle_col = 0,
  border_color = NA,
  fontsize_row = 6
)

# Select genes with opposite values in the two columns
opposite_genes <- merged_table[(merged_table$log2FoldChange_IL5 > 0 & merged_table$log2FoldChange_MepoCil < 0) | 
                                 (merged_table$log2FoldChange_IL5 < 0 & merged_table$log2FoldChange_MepoCil > 0), ]

# Print the list of genes with opposite values
print(opposite_genes$Gene)




# FDR < 0.01
common_genes_001 <- data.frame(
  Gene = c("KRT8", "TRIM2", "GAPDH", "CXXC5", "MFHAS1", "THRB", "VAV3", "CFL1", "ARHGAP21", 
           "NXN", "CCT5", "CHD7", "FAM172A", "STOML2", "CCDC28A", "UBB", "PCNX2", "CLPTM1L", 
           "GANC", "PDCL3", "DYRK1A", "MRPL47", "GPATCH8", "PSMD2", "WDR7", "BCL2L11", 
           "TEAD1", "ROBO1", "SH3RF1", "RAN", "TAB2", "PPID", "PCBD2", "GDA", "FAR2", 
           "MAP4K3", "MRPS6", "HTATIP2", "THADA", "FTH1", "DIAPH2", "MAGI3", "CDYL", 
           "PGAM1", "DUS1L", "EYA1", "NR3C1", "PHLPP1", "AKR1A1", "PARVA", "GLB1L", 
           "USP24", "SF1", "MYO5A", "ANOS1", "PDLIM1", "SELENOW", "ST8SIA4", "SMG7", 
           "MBTD1", "AAMP", "LDLR", "RPRD1A", "SULT1A1", "MID1IP1", "MBD5", "VAMP3", 
           "CLDN16", "ACVR2A", "STAG1", "HNF4G", "ATG10", "STK3", "CCDC89")
)

# Merge the tables to include only the genes in common
merged_table <- merge(common_genes_001, IL5t, by = "Gene")  # Add FoldChange1
merged_table <- merge(merged_table, cilDE, by = "Gene")  # Add FoldChange2


library(pheatmap)
# Prepare the data for the heat map
# Assuming `merged_table` contains columns: "Gene", "log2FoldChange_IL5", "log2FoldChange_MepoCil"
# Set row names as genes and select the two fold change columns
heatmap_data <- merged_table[, c("log2FoldChange_IL5", "log2FoldChange_MepoCil")]
rownames(heatmap_data) <- merged_table$Gene

# Generate the heat map
pheatmap(
  mat = heatmap_data,               # Data matrix for heat map
  scale = "none",                   # No scaling applied (use "row" or "column" for scaling)
  cluster_rows = TRUE,              # Cluster rows (genes)
  cluster_cols = TRUE,              # Cluster columns (conditions)
  color = colorRampPalette(c("royalblue", "white", "brown3"))(50), # Color gradient
  main = "DE IL-5 stim vs. Mepo (Ciliated cells), FDR<0.01", # Title of the heat map
  angle_col = 0,
  border_color = NA,
  fontsize_row = 6)

# Select genes with opposite values in the two columns
opposite_genes_001 <- merged_table[(merged_table$log2FoldChange_IL5 > 0 & merged_table$log2FoldChange_MepoCil < 0) | 
                                 (merged_table$log2FoldChange_IL5 < 0 & merged_table$log2FoldChange_MepoCil > 0), ]

# Print the list of genes with opposite values
print(opposite_genes_001$Gene)




# Venn diagrams
library(venneuler)

venn <- venneuler(c("IL-5"=5075, "Mepo/Cil"=348, "Mepo/Cil&IL-5"=111))
plot(venn, col = c("grey", "skyblue"))
title("IL-5 (FDA<0.05) vs Mepo, Ciliated cells")

venn <- venneuler(c("Mepo/Cil"=348, "IL-5"=5075, "Mepo/Cil&IL-5"=111))
plot(venn, col = c("grey", "cornflowerblue"))
title("IL-5 (FDA<0.05) vs Mepo, Ciliated cells")

venn2 <- venneuler(c("Mepo/Cil"=274, "IL-5"=3175, "Mepo/Cil&IL-5"=74))
plot(venn2, col = c("grey", "deeppink"))
title("IL-5 (FDA<0.01) vs Mepo, Ciliated cells")

venn3 <- venneuler(c("Epith."=1326, "Mepo/Cil"=348, "Mepo/Cil&IL-5"=21))
plot(venn3, col = c("red4", "skyblue"))
title("Epithelial vs Mepo, Ciliated cells")




# Calculate expression levels of the gene X per patient per condition
# Specify the gene of interest
gene_of_interest <- "IL5RA"

# Extract the RNA expression matrix
expression_matrix <- GetAssayData(merged_seurat.harmony2, assay = "RNA", slot = "data")

# Check if the gene is present in the expression matrix
if (!(gene_of_interest %in% rownames(expression_matrix))) {
  stop("Gene not found in the expression matrix.")
}

# Extract metadata containing 'Patient' and 'Type' information
metadata <- merged_seurat.harmony2@meta.data

# Verify the required columns exist in metadata
if (!all(c("Patient", "Type") %in% colnames(metadata))) {
  stop("'Patient' and/or 'Type' column not found in metadata.")
}

# Create a data frame with expression values, patient, and type
expression_data <- data.frame(
  Patient = metadata$Patient,
  Type = metadata$Type,
  Expression = expression_matrix[gene_of_interest, ]
)

# Calculate the average expression for IL5RA grouped by Patient and Type
average_expression <- aggregate(Expression ~ Patient + Type, data = expression_data, FUN = mean)

# Rename columns for clarity
colnames(average_expression) <- c("Patient", "Type", "Average_Expression")

# Print the resulting table
print(average_expression)

# Optionally save the table to a CSV file
# write.csv(average_expression, "IL5RA_average_expression_per_patient_type.csv", row.names = FALSE)

# Load ggplot2
library(ggplot2)

# Plot the average expression per patient for each Type
ggplot(average_expression, aes(x = Patient, y = Average_Expression, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +  # Grouped bar plot
  labs(
    title = "Average IL5RA Expression per Patient per Condition",
    x = "Patient",
    y = "Average Expression"
  ) +
  scale_fill_manual(values = c("A" = "burlywood", "B" = "aquamarine4")) +  # Customize colors for Type
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Centered and bold title
    axis.text.x = element_text(angle = 45, hjust = 1)      # Rotate x-axis labels for readability
  )





# Compare DE genes IL-5 stimulated and Basal Mepo treated 

#Load tables
basDE <- read.csv("~/Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/Res-200-5000-paired/basal30h.csv", header = T)
IL5t <- read.csv("~/Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/IL_5_resp_significant.csv", header = T) # From https://doi.org/10.1111/all.14297

# Modify the tables to get just two columns: Gene and log2FoldChange_IL5 or log2FoldChange_MepoCil
names(basDE)[1] <- "Gene" # Name the first column
names(basDE)[3] <- "log2FoldChange_MepoBas" # Rename the column Foldchange
basDE <- basDE[,-c(2,4,5,6,7,8)] # Remove other columns
names(IL5t)[2] <- "log2FoldChange_IL5"
IL5t <- IL5t[, -c(1,3,4,5)]

# This is a list of genes found to be DE in the table IL5 stimulation and in the Basal from IRIS study
common_genes <- data.frame(
  Gene = c("DHX32", "DPYSL3", "GDA", "HLA-A", "HLA-DRA", "PFKP", "RFX2", "SAMHD1", "SLC12A2")
)

# Merge the tables to include only the genes in common
merged_table <- merge(common_genes, IL5t, by = "Gene")  # Add FoldChange1
merged_table <- merge(merged_table, basDE, by = "Gene")  # Add FoldChange2


library(pheatmap)
# Prepare the data for the heat map
# Assuming `merged_table` contains columns: "Gene", "log2FoldChange_IL5", "log2FoldChange_MepoBas"
# Set row names as genes and select the two fold change columns
heatmap_data <- merged_table[, c("log2FoldChange_IL5", "log2FoldChange_MepoBas")]
rownames(heatmap_data) <- merged_table$Gene

# Generate the heat map
pheatmap(
  mat = heatmap_data,               # Data matrix for heat map
  scale = "none",                   # No scaling applied (use "row" or "column" for scaling)
  cluster_rows = TRUE,              # Cluster rows (genes)
  cluster_cols = TRUE,              # Cluster columns (conditions)
  color = colorRampPalette(c("royalblue", "white", "brown3"))(50), # Color gradient
  main = "DE IL-5 stim vs. Mepo (Basal cells)", # Title of the heat map
  angle_col = 0,
  border_color = NA,
  fontsize_row = 6)




# Compare DE genes IL-5 stimulated and Ionocytes Mepo treated 

#Load tables
ionoDE <- read.csv("~/Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/Res-200-5000-paired/iono10h.csv", header = T)
IL5t <- read.csv("~/Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/IL_5_resp_significant.csv", header = T) # From https://doi.org/10.1111/all.14297

# Modify the tables to get just two columns: Gene and log2FoldChange_IL5 or log2FoldChange_MepoCil
names(ionoDE)[1] <- "Gene" # Name the first column
names(ionoDE)[3] <- "log2FoldChange_MepoIono" # Rename the column Foldchange
ionoDE <- ionoDE[,-c(2,4,5,6,7,8)] # Remove other columns
names(IL5t)[2] <- "log2FoldChange_IL5"
IL5t <- IL5t[, -c(1,3,4,5)]

# This is a list of genes found to be DE in the table IL5 stimulation and in the Ionocytes  from IRIS study
common_genes <- data.frame(
  Gene =c("ATP1A1", "CFL1", "GAPDH", "HLA-A", "HLA-B", "KRT7", "KRT8", "SEC11C", "SPINT2")
)

# Merge the tables to include only the genes in common
merged_table <- merge(common_genes, IL5t, by = "Gene")  # Add FoldChange1
merged_table <- merge(merged_table, ionoDE, by = "Gene")  # Add FoldChange2


library(pheatmap)
# Prepare the data for the heat map
# Assuming `merged_table` contains columns: "Gene", "log2FoldChange_IL5", "log2FoldChange_MepoIono"
# Set row names as genes and select the two fold change columns
heatmap_data <- merged_table[, c("log2FoldChange_IL5", "log2FoldChange_MepoIono")]
rownames(heatmap_data) <- merged_table$Gene

# Generate the heat map
pheatmap(
  mat = heatmap_data,               # Data matrix for heat map
  scale = "none",                   # No scaling applied (use "row" or "column" for scaling)
  cluster_rows = TRUE,              # Cluster rows (genes)
  cluster_cols = TRUE,              # Cluster columns (conditions)
  color = colorRampPalette(c("royalblue", "white", "brown3"))(50), # Color gradient
  main = "DE IL-5 stim vs. Mepo (Ionocytes)", # Title of the heat map
  angle_col = 0,
  border_color = NA,
  fontsize_row = 6)



# Evaluate the log2FC of the genes DE in ciliated and SMAD2 dependent - SMAD2 goes down in Cil cells post MEPO
DE_cil<-read.csv("~/Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/Res-200-5000-paired/ciliated50h.csv")
names(DE_cil)[1] <- "Gene"
gene_list <- c("ADK", "ALCAM", "ARHGAP21", "BCL2L11", "BMP2K", "CDYL", "CREB3", "CXXC5",
               "DENND1B", "DMTF1", "DYRK1A", "GADD45A", "GANC", "GUSBP1", "KRT8", "MAGI3",
               "MAP3K4", "MAPK8", "MBTD1", "PDLIM1", "PHF21A", "PHLPP1", "PKM", "PLAG1",
               "POGZ", "RAB4A", "RARB", "RASA2", "RBMS3", "SETBP1", "SH3RF1", "SMG7",
               "TALDO1", "THADA", "TMEM87A", "TRIM8", "VPS13B")

# Load required libraries
library(ggplot2)

# Filter DE_cil for genes present in gene_list
DE_cil_subset <- DE_cil[DE_cil$Gene %in% gene_list, ]

# Check if any genes matched the filter
if (nrow(DE_cil_subset) == 0) {
  stop("No matching genes from gene_list found in DE_cil.")
}

# Order genes by log2FoldChange from lowest to highest
DE_cil_subset <- DE_cil_subset[order(DE_cil_subset$log2FoldChange), ]

# Create a bar plot with flipped axes and ordered genes
ggplot(DE_cil_subset, aes(x = reorder(Gene, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("TRUE" = "brown3", "FALSE" = "royalblue3"),
                    name = "Direction",
                    labels = c("Downregulated", "Upregulated")) +
  coord_flip() +  # Flip the axes
  theme_minimal() +
  labs(title = "Log2 Fold Change SMAD2-dependent genes, DE Ciliated, post-Mepo",
       x = "Gene",
       y = "Log2 Fold Change") +
  theme(axis.text.y = element_text(size = 10))  # Adjust text size for better readability


# Here is the list of DEGs from Ciliated cells, post-Mepo that are listed to be Cell-Cell Junction related
DE_cil<-read.csv("~/Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/Res-200-5000-paired/ciliated50h.csv")
names(DE_cil)[1] <- "Gene"
overlapping_genes <- c("CFL1","CLDN16","KRT8","MAGI3","PRKD1")

# Filter DE_cil for genes present in gene_list
DE_cil_subset <- DE_cil[DE_cil$Gene %in% overlapping_genes, ]

# Check if any genes matched the filter
if (nrow(DE_cil_subset) == 0) {
  stop("No matching genes from gene_list found in DE_cil.")
}

# Order genes by log2FoldChange from lowest to highest
DE_cil_subset <- DE_cil_subset[order(DE_cil_subset$log2FoldChange), ]

# Create a bar plot with flipped axes and ordered genes
ggplot(DE_cil_subset, aes(x = reorder(Gene, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("TRUE" = "brown3", "FALSE" = "royalblue3"),
                    name = "Direction",
                    labels = c("Downregulated", "Upregulated")) +
  coord_flip() +  # Flip the axes
  theme_minimal() +
  labs(title = "Log2 Fold Change Cell-Cell Junction, DE Ciliated, post-Mepo",
       x = "Gene",
       y = "Log2 Fold Change") +
  theme(axis.text.y = element_text(size = 10))  # Adjust text size for better readability


# Remodeling genes
# Here is the list of DEGs from Ciliated cells, post-Mepo that are listed to be Remodelling related

remodelling_genes <- c("ABI1", "ARHGAP21", "BSG", "CFL1", "CHD7", "CNN3", "GAPDH", "GSK3B",
                       "MAP2K4", "MDK", "NR3C1", "NR3C2", "PLTP", "PRKD1", "RARB", "SMAD2",
                       "SMARCAD1", "STK3", "TIMP1", "TLK1"
)


# Filter DE_cil for genes present in gene_list
DE_cil_subset <- DE_cil[DE_cil$Gene %in% remodelling_genes, ]

# Check if any genes matched the filter
if (nrow(DE_cil_subset) == 0) {
  stop("No matching genes from gene_list found in DE_cil.")
}

# Order genes by log2FoldChange from lowest to highest
DE_cil_subset <- DE_cil_subset[order(DE_cil_subset$log2FoldChange), ]

# Create a bar plot with flipped axes and ordered genes
ggplot(DE_cil_subset, aes(x = reorder(Gene, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("TRUE" = "brown3", "FALSE" = "royalblue3"),
                    name = "Direction",
                    labels = c("Downregulated", "Upregulated")) +
  coord_flip() +  # Flip the axes
  theme_minimal() +
  labs(title = "Log2 Fold Change Remodelling genes, DE Ciliated, post-Mepo",
       x = "Gene",
       y = "Log2 Fold Change") +
  theme(axis.text.y = element_text(size = 10))  # Adjust text size for better readability

# DE SAMD2-dependend found in the IL-5 signature, cil cells, post-mepo list of genes with opposite change compared to IL-5 stimulation

signature_genes_SMAD2dep <- c("BCL2L11","CXXC5","DENND1B","DMTF1","GANC","KRT8","MAGI3","MBTD1","PDLIM1")

# Filter DE_cil for genes present in gene_list
DE_cil_subset <- DE_cil[DE_cil$Gene %in% signature_genes_SMAD2dep, ]

# Check if any genes matched the filter
if (nrow(DE_cil_subset) == 0) {
  stop("No matching genes from gene_list found in DE_cil.")
}

# Order genes by log2FoldChange from lowest to highest
DE_cil_subset <- DE_cil_subset[order(DE_cil_subset$log2FoldChange), ]

# Create a bar plot with flipped axes and ordered genes
ggplot(DE_cil_subset, aes(x = reorder(Gene, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("TRUE" = "brown3", "FALSE" = "royalblue3"),
                    name = "Direction",
                    labels = c("Downregulated", "Upregulated")) +
  coord_flip() +  # Flip the axes
  theme_minimal() +
  labs(title = "Log2 Fold Change SMAD2-dependent-IL5-signature",
       x = "Gene",
       y = "Log2 Fold Change") +
  theme(axis.text.y = element_text(size = 10))  # Adjust text size for better readability

# Remodeling genes
# Here is the list of DEGs from Ciliated cells, post-Mepo that are listed to be remodeling related

remodelling_genes <- c("ABI1", "ARHGAP21", "CDK2AP1", "CFL1", "CHD7", "CNN3",
  "GAPDH", "GSK3B", "MAP2K4", "MDK", "NR3C1", "NR3C2",
  "PLTP", "PRKD1", "PSMD2", "RARB", "SMAD2", "SMARCAD1",
  "STK3", "TIMP1", "TLK1", "VIM"
)


# Filter DE_cil for genes present in gene_list
DE_cil_subset <- DE_cil[DE_cil$Gene %in% remodelling_genes, ]

# Check if any genes matched the filter
if (nrow(DE_cil_subset) == 0) {
  stop("No matching genes from gene_list found in DE_cil.")
}

# Order genes by log2FoldChange from lowest to highest
DE_cil_subset <- DE_cil_subset[order(DE_cil_subset$log2FoldChange), ]

# Create a bar plot with flipped axes and ordered genes
ggplot(DE_cil_subset, aes(x = reorder(Gene, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("TRUE" = "brown3", "FALSE" = "royalblue3"),
                    name = "Direction",
                    labels = c("Downregulated", "Upregulated")) +
  coord_flip() +  # Flip the axes
  theme_minimal() +
  labs(title = "Log2 Fold Change Remodeling genes, DE Ciliated, post-Mepo",
       x = "Gene",
       y = "Log2 Fold Change") +
  theme(axis.text.y = element_text(size = 10))  # Adjust text size for better readability

# Evaluate the log2FC of the genes DE in ciliated and related to Rhinovirus infection - SMAD2 goes down in Cil cells post MEPO
DE_cil<-read.csv("~/Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/Res-200-5000-paired/res_ciliated_keep50.csv")
names(DE_cil)[1] <- "Gene"
gene_list <- c("HLA-C","LDLR", "MAPK8", "MYO5A", "NR3C1")

# Load required libraries
library(ggplot2)

# Filter DE_cil for genes present in gene_list
DE_cil_subset <- DE_cil[DE_cil$Gene %in% gene_list, ]

# Check if any genes matched the filter
if (nrow(DE_cil_subset) == 0) {
  stop("No matching genes from gene_list found in DE_cil.")
}

# Order genes by log2FoldChange from lowest to highest
DE_cil_subset <- DE_cil_subset[order(DE_cil_subset$log2FoldChange), ]

# Create a bar plot with flipped axes and ordered genes
ggplot(DE_cil_subset, aes(x = reorder(Gene, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("TRUE" = "brown3", "FALSE" = "royalblue3"),
                    name = "Direction",
                    labels = c("Downregulated", "Upregulated")) +
  coord_flip() +  # Flip the axes
  theme_minimal() +
  labs(title = "Log2 Fold Change Rhinovirus Infection genes, DE Ciliated, post-Mepo",
       x = "Gene",
       y = "Log2 Fold Change") +
  theme(axis.text.y = element_text(size = 10))  # Adjust text size for better readability


# IL-13 Check #

#Load tables
cilDE <- read.csv("~/Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/Res-200-5000-paired/ciliated50h.csv", header = T)
IL13t <- read.csv("~/Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/IL13_signature_significant_ciliated.csv", header = T) # From https://doi.org/10.1111/all.14297

# Modify the tables to get just two columns: Gene and log2FoldChange_IL5 or log2FoldChange_MepoCil
names(cilDE)[1] <- "Gene" # Name the first column
names(cilDE)[3] <- "log2FoldChange_MepoCil" # Rename the column Foldchange
cilDE <- cilDE[,-c(2,4,5,6,7,8)] # Remove other columns
IL13t <- IL13t[, -c(7,3,4,5,6)] # Remove other columns
names(IL13t)[2] <- "log2FoldChange_IL13"

# This is a list of genes found to be DE in the table IL13 stimulation and in the ciliated cells from IRIS study
common_genes <- data.frame(
  Gene = genes <- c("AKR1A1",
                    "C3orf52",
                    "GAPDH",
                    "MAGI3",
                    "RASA2",
                    "SRBD1",
                    "TIMP1"
  )
  
)


merged_table <- merge(common_genes, IL13t, by = "Gene")  # Add FoldChange1
merged_table <- merge(merged_table, cilDE, by = "Gene")  # Add FoldChange2


library(pheatmap)
# Prepare the data for the heat map
# Assuming `merged_table` contains columns: "Gene", "log2FoldChange_IL13", "log2FoldChange_MepoCil"
# Set row names as genes and select the two fold change columns
heatmap_data <- merged_table[, c("log2FoldChange_IL13", "log2FoldChange_MepoCil")]
rownames(heatmap_data) <- merged_table$Gene

# Generate the heat map
# Determine the limits for the color scale
heatmap_range <- max(abs(heatmap_data), na.rm = TRUE)

# pheatmap(
#   mat = heatmap_data,
#   scale = "none",
#   cluster_rows = TRUE,
#   cluster_cols = TRUE,
#   color = colorRampPalette(c("blue", "white", "red"))(100),
#   breaks = seq(-heatmap_range, heatmap_range, length.out = 101),
#   main = "DE IL-5 stim vs. Mepo (Ciliated cells)",
#   angle_col = 0,
#   border_color = NA,
#   fontsize_row = 6
# )
# Cap the heatmap data at -1 and 1
heatmap_data_capped <- heatmap_data
heatmap_data_capped[heatmap_data_capped < -1] <- -1
heatmap_data_capped[heatmap_data_capped > 1] <- 1

# Define color and breaks
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- seq(-1, 1, length.out = 101)

# Draw the heatmap
pheatmap(
  mat = heatmap_data_capped,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = color_palette,
  breaks = breaks,
  main = "DE IL-13 stim vs. Mepo (Ciliated cells)",
  angle_col = 0,
  border_color = NA,
  fontsize_row = 6
)
