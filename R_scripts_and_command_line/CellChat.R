# CellCHat
c("#606372FF", "#79A8A4FF", "#B2AD8FFF", "steelblue", "#DEC18CFF", "#92A185FF")

library(CellChat)
library(Seurat)
library(SeuratObject)
library(dplyr)

# Create a CellChat object from Seurat object
merged_seurat.harmony <- readRDS("~/Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/merged_seurat.harmonyCORRECTED.rds")

# Create a CellChat object
merged_seurat.harmony@meta.data$samples <- merged_seurat.harmony@meta.data$Type
cellchat <- createCellChat(object = merged_seurat.harmony, group.by = "annotation1")

# Set the ligand-receptor interaction database (default: for human)
CellChatDB <- CellChatDB.human # Use CellChatDB.mouse for mouse datasets
cellchat@DB <- CellChatDB

# View the CellChat object
cellchat

# Subset relevant ligand-receptor interactions
cellchat <- subsetData(cellchat) 

# Identify overexpressed genes and interactions
# devtools::install_github('immunogenomics/presto')
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Compute communication probabilities
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10) # Adjust min.cells as needed

# Infer communication pathways
cellchat <- computeCommunProbPathway(cellchat)

# Aggregate communication networks
cellchat <- aggregateNet(cellchat)

# View summary of interactions
print(cellchat)

# Network plot
par(mfrow = c(1,2), xpd = TRUE )
netVisual_circle(cellchat@net$count, vertex.weight = TRUE, weight.scale = TRUE, label.edge= FALSE, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = TRUE, weight.scale = TRUE, label.edge= FALSE, title.name = "Interaction weights/strength")

dev.off() # To reset the setting of the multiple plots on the same area

# Visualisation of the pathways involved

cellchat@netP[["pathways"]] # It shows the list of the pathways involved into the cell-cell communication network

extractEnrichedLR(cellchat, signaling = c(cellchat@netP[["pathways"]]),  
                  geneLR.return = TRUE) # This returns the list of genes involved


pathways.show <- c("COLLAGEN") 
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle") # layout "hierarchy", "circle", "chord" or "spatial"

# Contribution of each signaling pathway
netAnalysis_contribution(cellchat, signaling = c(cellchat@netP[["pathways"]]),
                         title = "Contribution of each LR")

netAnalysis_contribution(cellchat, signaling = c(cellchat@netP[["pathways"]][1:10]),
                         title = "Contribution top 10 LR")

extractEnrichedLR(cellchat, signaling = "JAM", geneLR.return = FALSE)
netAnalysis_contribution(cellchat, signaling = "JAM")

# Heatmap of communication strength
netVisual_heatmap(cellchat)

# Signaling role analysis
cellchat <- netAnalysis_computeCentrality(cellchat)
netAnalysis_computeCentrality(cellchat)
netAnalysis_signalingRole_scatter(cellchat)


# Save the first, common cellchat object
saveRDS(cellchat, "~/Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/cellchat_general_CORRECTED.rds")
################################################################################

# Compare A and B
library(CellChat)
library(patchwork)
library(NMF)
library(ggalluvial)
options(stringsAsFactors = FALSE)

merged_seurat.harmony
table(merged_seurat.harmony[["Type"]])
data.input <- GetAssayData(merged_seurat.harmony, layer = "data")
meta <- merged_seurat.harmony@meta.data

cell.use <- rownames(meta)[meta$samples == "A"]
data.input <- data.input[, cell.use]
meta <- meta[meta$samples == "A", ]

cellchatA <- createCellChat(object = data.input, meta = meta, group.by = "annotation1")

# Set the ligand-receptor interaction database (default: for human)
CellChatDB <- CellChatDB.human # Use CellChatDB.mouse for mouse datasets
cellchatA@DB <- CellChatDB

# A
# Subset relevant ligand-receptor interactions
cellchatA <- subsetData(cellchatA) 

# Identify overexpressed genes and interactions
# devtools::install_github('immunogenomics/presto')
cellchatA <- identifyOverExpressedGenes(cellchatA)
cellchatA <- identifyOverExpressedInteractions(cellchatA)

# Compute communication probabilities
cellchatA <- computeCommunProb(cellchatA)
cellchatA <- filterCommunication(cellchatA, min.cells = 10) # Adjust min.cells as needed

cellchatA <- computeCommunProbPathway(cellchatA)

cellchatA <- aggregateNet(cellchatA)

cellchatA <- netAnalysis_computeCentrality(cellchatA, slot.name = "netP")

# Identify and visualize outgoing and incoming communication patterns of target cells
selectK(cellchatA, pattern = "outgoing")
nPatterns = 2# This number comes from the plot generated from the previous command. Choose the number after which you see the first drop
dev.off()
chellchatA <- identifyCommunicationPatterns(cellchatA, pattern = "outgoing", k = nPatterns, width = 5, height = 9)

selectK(cellchatA, pattern = "incoming")
nPatterns = 2 # This number comes from the plot generated from the previous command. Choose the number after which you see the first drop
dev.off()
chellchatA <- identifyCommunicationPatterns(cellchatA, pattern = "incoming", k = nPatterns, width = 5, height = 9)

# Identify signaling groups based on functional similarity
cellchatA <- computeNetSimilarity(cellchatA, type = "functional")
# Before going further, you need to force a different version of python
py_install("umap-learn")
py_install("umap-learn", envname = "r-reticulate-py310")

# ! Restart R first: Session > Restart R
library(reticulate)
use_condaenv("r-reticulate-py311", required = TRUE)
py_config()  # Verify that Python 3.11 is being used

py_module_available("umap") # If you get TRUE as output, then you can go on

cellchatA <- netEmbedding(cellchatA, type = "functional")
cellchatA <- netClustering(cellchatA, type = "functional", do.parallel = FALSE)

# Identify signaling groups based on structure similarity
cellchatA <- computeNetSimilarity(cellchatA, type = "structural")
cellchatA <- netEmbedding(cellchatA, type = "structural")
cellchatA <- netClustering(cellchatA, type = "structural", do.parallel = FALSE)

# Save cellchat object for A
saveRDS(cellchatA, "~/Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/cellchat_analysed_A_CORRECTED.rds")

# B
table(merged_seurat.harmony[["Type"]])
data.input <- GetAssayData(merged_seurat.harmony, layer = "data")
meta <- merged_seurat.harmony@meta.data

cell.use <- rownames(meta)[meta$samples == "B"]
data.input <- data.input[, cell.use]
meta <- meta[meta$samples == "B", ]

cellchatB <- createCellChat(object = data.input, meta = meta, group.by = "annotation1")


# Set the ligand-receptor interaction database (default: for human)
CellChatDB <- CellChatDB.human 
cellchatB@DB <- CellChatDB
# Subset relevant ligand-receptor interactions
cellchatB <- subsetData(cellchatB) 

# Identify overexpressed genes and interactions
# devtools::install_github('immunogenomics/presto')
cellchatB <- identifyOverExpressedGenes(cellchatB)
cellchatB <- identifyOverExpressedInteractions(cellchatB)

# Compute communication probabilities
cellchatB <- computeCommunProb(cellchatB)
cellchatB <- filterCommunication(cellchatB, min.cells = 10) # Adjust min.cells as needed

cellchatB <- computeCommunProbPathway(cellchatB)

cellchatB <- aggregateNet(cellchatB)

cellchatB <- netAnalysis_computeCentrality(cellchatB, slot.name = "netP")

# Identify and visualize outgoing and incoming communication patterns of target cells
selectK(cellchatB, pattern = "outgoing")
nPatterns = 2 # This number comes from the plot generated from the previous command. Choose the number after which you see the first drop
chellchatB <- identifyCommunicationPatterns(cellchatB, pattern = "outgoing", k = nPatterns, width = 5, height = 9)

selectK(cellchatB, pattern = "incoming")
nPatterns = 3 # This number comes from the plot generated from the previous command. Choose the number after which you see the first drop
chellchatB <- identifyCommunicationPatterns(cellchatB, pattern = "incoming", k = nPatterns, width = 5, height = 9)

# Identify signaling groups based on functional similarity
cellchatB <- computeNetSimilarity(cellchatB, type = "functional")
cellchatB <- netEmbedding(cellchatB, type = "functional")
cellchatB <- netClustering(cellchatB, type = "functional", do.parallel = FALSE)

# Identify signaling groups based on structure similarity
cellchatB <- computeNetSimilarity(cellchatB, type = "structural")
cellchatB <- netEmbedding(cellchatB, type = "structural")
cellchatB <- netClustering(cellchatB, type = "structural", do.parallel = FALSE)

# Save cellchat object for B
saveRDS(cellchatB, "~/Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/cellchat_analysed_B_CORRECTED.rds")

# Load the cellchat objects
library(CellChat)
library(Seurat)
library(SeuratObject)
library(dplyr)
library(patchwork)
library(NMF)
library(ggalluvial)
options(stringsAsFactors = FALSE)

cellchatA <- readRDS("~/Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/cellchat_analysed_A_CORRECTED.rds")
cellchatB <- readRDS("~/Desktop/scRNA-seq_Brush1/IRIS/Full_dataset/cellchat_analysed_B_CORRECTED.rds")

object.list <- list(pre = cellchatA, post = cellchatB)
names(object.list)

cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat

# Compare difference s between A and B

# Compare the overall information flow of each signaling pathway
rankNet(cellchat, mode = "comparison", stacked = TRUE, do.stat = TRUE, color.use = c("grey80", "steelblue")) #"" "skyblue4"
rankNet(cellchat, mode = "comparison", stacked = FALSE, do.stat = TRUE, color.use = c("grey80", "steelblue"))

# Compare the total number of interactions and interaction strength
compareInteractions(cellchat, show.legend =  FALSE, group = c(1, 2), measure = "count", color.use = c("grey80", "steelblue"))
compareInteractions(cellchat, show.legend =  FALSE, group = c(1, 2), measure = "weight", color.use = c("grey80", "steelblue"))

# List the pathways involved per condition
cellchat@netP[["pre"]][["pathways"]]
cellchat@netP[["post"]][["pathways"]]

# List the genes related to the pathways involved per condition
extractEnrichedLR(cellchat, signaling = c(cellchat@netP[["pre"]][["pathways"]]),  
                                          geneLR.return = TRUE)

extractEnrichedLR(cellchat, signaling = c(cellchat@netP[["post"]][["pathways"]]),  
                  geneLR.return = TRUE)

## Scatterplot
# Compare outgoing/incoming interaction strength fall all the cell types
count.sum <- sapply(object.list, function(x) {
  rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)
})

# This controls the scaling of symbol sizes
weight.MinMax <- c(min(count.sum), max(count.sum))

# Determine axis ranges across all objects
all.x <- all.y <- c()
for (x in object.list) {
  df <- netAnalysis_signalingRole_scatter(x, weight.MinMax = weight.MinMax)$data
  all.x <- c(all.x, df$x)
  all.y <- c(all.y, df$y)
}
axis.lims <- list(x = range(all.x, na.rm = TRUE),
                  y = range(all.y, na.rm = TRUE))

# Generate plots with fixed scales
gg <- list()
for (i in seq_along(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]],
                                               title = names(object.list)[i],
                                               weight.MinMax = weight.MinMax
  ) +
    xlim(axis.lims$x) +
    ylim(axis.lims$y)
}

patchwork::wrap_plots(plots = gg)

# To limit the scatter plots to a specific pathway
pathway_of_interest <- "FN1"  # change to any pathway name, e.g. "JAM", "FN1", "POSTN"

# List pathways per condition (no change needed here — just for reference)
cellchat@netP[["pre"]][["pathways"]]
cellchat@netP[["post"]][["pathways"]]

# Filter enriched LR pairs to your pathway only
extractEnrichedLR(cellchat, signaling = pathway_of_interest, geneLR.return = TRUE)

# Scatterplot — filter count matrix to your pathway only before summing
count.sum <- sapply(object.list, function(x) {
  # Extract the count matrix for your pathway only
  net <- x@netP$prob[, , pathway_of_interest]  # 3D array: sender x receiver x pathway
  rowSums(net) + colSums(net) - diag(net)
})

weight.MinMax <- c(min(count.sum), max(count.sum))

all.x <- all.y <- c()
for (x in object.list) {
  df <- netAnalysis_signalingRole_scatter(
    x,
    signaling = pathway_of_interest,   # <-- key argument added here
    weight.MinMax = weight.MinMax
  )$data
  all.x <- c(all.x, df$x)
  all.y <- c(all.y, df$y)
}

axis.lims <- list(x = range(all.x, na.rm = TRUE),
                  y = range(all.y, na.rm = TRUE))

gg <- list()
for (i in seq_along(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(
    object.list[[i]],
    signaling   = pathway_of_interest,  # <-- key argument added here
    title       = names(object.list)[i],
    weight.MinMax = weight.MinMax
  ) +
    xlim(axis.lims$x) +
    ylim(axis.lims$y)
}

patchwork::wrap_plots(plots = gg)


# Identify signalling changes associated with one cell type
netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Basal", color.use = c("black", "grey80", "steelblue"))
netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Ciliated", color.use = c("black", "grey80", "steelblue"))
netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Deuterosomal", color.use = c("black", "grey80", "steelblue"))
netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Mast cells", color.use = c("black", "grey80", "steelblue"))
netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Ionocytes", color.use = c("black", "grey80", "steelblue"))
netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Goblet", color.use = c("black", "grey80", "steelblue"))
netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Club", color.use = c("black", "grey80", "steelblue"))
netAnalysis_signalingChanges_scatter(cellchat, idents.use = "T-cells", color.use = c("black", "grey80", "steelblue"))
netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Macrophages", color.use = c("black", "grey80", "steelblue"))
netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Dendritic", color.use = c("black", "grey80", "steelblue"))
netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Mucociliated", color.use = c("black", "grey80", "steelblue"))


# Print specific interaction of a certain pathway
extractEnrichedLR(cellchat, signaling = "LAMININ", geneLR.return = T)

## Circle plots
# Show the number of interactions between any two cell populations
# compute the maximum number of cells and the maximum number of interactions
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))

library(CellChat)

par(mfrow = c(1, 2), xpd = TRUE)

for (i in seq_along(object.list)) {
  # Extract unique cell identities
  cell_types <- unique(object.list[[i]]@idents)
  
  # Generate colors using CellChat's default scPalette
  color.use <- setNames(scPalette(length(cell_types)), cell_types)
  
  # Debugging: Print to verify correct color mapping
  print(color.use)
  
  # Generate the circle plot
  netVisual_circle(
    object.list[[i]]@net$count,
    weight.scale = TRUE,
    label.edge = FALSE,
    edge.weight.max = weight.max[2],
    edge.width.max = 12,
    arrow.size = 0.05,
    color.use = color.use,  # Using CellChat colors
    title.name = paste0("Number of interactions - ", names(object.list)[i])
  )
}

# For a selected pathway
weight.max <- getMaxWeight(object.list, slot.name = "netP", attribute = "JAM")

# Ensure plots do not overlap
par(mfrow = c(1, 2), mar = c(4, 4, 2, 2), xpd = TRUE)

for (i in seq_along(object.list)) {  
  netVisual_aggregate(
    object.list[[i]], 
    signaling = "JAM", 
    layout = "circle",
    edge.weight.max = weight.max[1], 
    edge.width.max = 10, 
    arrow.size = 0.05,
    signaling.name = paste("JAM -", names(object.list)[i])  # Corrected title
  )
}

# Reset plotting parameters (optional)
par(mfrow = c(1, 1))

# Show differential number of interactions or interaction strength among
# different cell populations, red (increased signaling), blue (decreased signaling)
par(mfrow = c(1,2), xpd = TRUE)
netVisual_diffInteraction(cellchat, comparison = c(1,2), measure = "count",
                          weight.scale = TRUE, arrow.size = 0.1)
netVisual_diffInteraction(cellchat, comparison = c(1,2), measure = "weight",
                          weight.scale = TRUE, arrow.size = 0.1)

# Simplify to the cell type level
group.cellType <- c(rep("Basal", 3), rep("Ciliated", 3), rep("Ionocytes", 3))
group.cellType <- factor(group.cellType, levels = c("Basal", "Ciliated", "Ionocytes"))
object.list <- lapply(object.list, function(x) {
                  mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"),
                           attribute = c("idents", "count", "count.merged"))

# show the number of interactions or interaction strength
par(mfrow = c(1,2), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = TRUE,
                   label.edge = TRUE, edge.weight.max = weight.max[3], edge.width.max = 12,
                   arrow.size = 0.5,
                   title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

par(mfrow = c(1,2), xpd = TRUE)
netVisual_diffInteraction(cellchat, weight.scale = TRUE, comparison = c(1,2),
                          arrow.size = 0.5, measure = "count.merged", label.edge = TRUE)
netVisual_diffInteraction(cellchat, weight.scale = TRUE, comparison = c(1,2),
                          arrow.size = 0.5, measure = "weight.merged", label.edge = TRUE)


## Heatmaps
all_pathways <- union(object.list[[1]]@netP$pathways,
                      object.list[[2]]@netP$pathways)

ht1 <- netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "all",
                                         signaling = all_pathways, title = names(object.list)[1],
                                         width =8, height = 25, color.heatmap = "PuBuGn")
ht2 <- netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "all",
                                         signaling = all_pathways, title = names(object.list)[2],
                                         width =8, height = 25, color.heatmap = "PuBuGn")
# draw(ht1, ht2, ht_gap = unit(0.5, "cm")) > Didn't work, so I used the following for plotting;
library(ComplexHeatmap)
ht_list <- ht1 + ht2
draw(ht_list, ht_gap = unit(0.5, "cm"))

# Plot separately outgoing and incoming signaling
ht_3 <- netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "outgoing",
                                          signaling = all_pathways, title = names(object.list)[1],
                                          width =5, height = 20, color.heatmap = "PuRd")
ht_4 <- netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "outgoing",
                                          signaling = all_pathways, title = names(object.list)[2],
                                          width =5, height = 20, color.heatmap = "PuRd")
ht_list <- ht_3 + ht_4
draw(ht_list, ht_gap = unit(0.5, "cm"))

ht_5 <- netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "incoming",
                                          signaling = all_pathways, title = names(object.list)[1],
                                          width =5, height = 20, color.heatmap = "BuGn")
ht_6 <- netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "incoming",
                                          signaling = all_pathways, title = names(object.list)[2],
                                          width =5, height = 20, color.heatmap = "BuGn")
ht_list <- ht_5 + ht_6
draw(ht_list, ht_gap = unit(0.5, "cm"))

# Selected pathways
par(mfrow = c(1,2), xpd = TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = c("LAMININ"),
                               title.name = paste("LAM", "signaling", names(object.list)[i]),
                               color.heatmap = "Reds")
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

# Show differential interaction number & interaction strength using heatmap
gg1 <- netVisual_heatmap(cellchat, comparison =  c(1,2), measure = "count")
gg2 <- netVisual_heatmap(cellchat, comparison =  c(1,2), measure = "weight")
gg1 + gg2

## Bubble plots
# compare communication probabilities mediated by ligand-recepor pairs
# all sources and all targets
netVisual_bubble(cellchat, comparison = c(1,2), angle.x = 45, color.text = c("grey80", "steelblue")) 
# from selected sources and targets cell groups
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1, 2, 7), # source.use = the cell type you want to compare, target.use = the cell types to compare to source.use
                  comparison = c(1,2), angle.x = 45, color.text = c("grey80", "steelblue")) 

# identify the up-regulated ligand-receptor pairs
netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1:11),
                 comparison = c(1,2), max.dataset = 2,
                 title.name = "Increased signaling post Mepo", angle.x = 45,
                 remove.isolate = TRUE, color.text = c("grey80", "steelblue"))

# identify the down-regulated ligand-receptor pairs
netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1:11),
                 comparison = c(1,2), max.dataset = 1,
                 title.name = "Decreased signaling post Mepo", angle.x = 45,
                 remove.isolate = TRUE, color.text = c("grey80", "steelblue"))
## Violin plot
View(cellchat@meta)

plotGeneExpression(cellchat, signaling = "FLRT", split.by = "samples",
                   color.use = c("grey80", "steelblue")) + theme(legend.position = "right")


# Calculate the aggregated cell-cell communication network.
# Due to the complicated cell-cell communication network, we can examine the signaling sent from each cell group. Here we also control the parameter edge.weight.max so that we can compare edge weights between differet networks.
# A
mat <- cellchatA@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
     mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
     mat2[i, ] <- mat[i, ]
     netVisual_circle(mat2, vertex.weight = TRUE, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i], arrow.width = 0.01, arrow.size = 0.003)
}
dev.off()
# B
mat <- cellchatB@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = TRUE, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i], arrow.width = 0.01, arrow.size = 0.003)
}
dev.off()


###
# Contribution top 10 pathways A vs B

netAnalysis_contribution(cellchatA, signaling = c(cellchatA@netP[["pathways"]][1:10]),
                         title = "Contribution top 10 LR, pre")
netAnalysis_contribution(cellchatB, signaling = c(cellchatB@netP[["pathways"]][1:10]),
                         title = "Contribution top 10 LR, post")



