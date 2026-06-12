suppressPackageStartupMessages({
    library(CellChat)
    library(RColorBrewer)
    library(ggplot2)
    library(dplyr)
    library(tidyverse)
    library(ggraph)
    library(wesanderson)
    library(ggrepel)
    library(Matrix)
})

# Load and catenate metadata
myeloid <- readRDS("../finalObjects_v3/myeloid_tissue_meta.Rds")  %>% select(cell, dataset, sample, Type, final_annotation) %>% mutate(label = "Myeloid_tissue")
b <- readRDS("../finalObjects_v3/b_tissue_meta.Rds") %>% select(cell, sample, dataset, Type, final_annotation) %>% mutate(label = "Bcell_tissue")
tnk <- readRDS("../finalObjects_v3/t_nk_tissue_meta.Rds")  %>% select(cell, sample, dataset, Type, final_annotation) %>% mutate(label = "TNK_tissue")
glom <- readRDS("../finalObjects_v3/glom_meta.Rds") %>% select(cell, sample, dataset, Type, final_annotation) %>% mutate(label = "GLOM")

all_meta <- bind_rows(
  myeloid, b, tnk, glom,
) %>% filter(dataset == "scRNAseq")

rm(myeloid, b, tnk, glom)
gc()

# Load and catenate normalized expression data
myeloidNorm <- readRDS("../finalObjects_v3/myeloid_tissue_norm.Rds")
myeloidNorm <- myeloidNorm[, colnames(myeloidNorm) %in% all_meta$cell[all_meta$label == "Myeloid_tissue"]]

bNorm <- readRDS("../finalObjects_v3/b_tissue_norm.Rds")
bNorm <- bNorm[, colnames(bNorm) %in% all_meta$cell[all_meta$label == "Bcell_tissue"]]

tnkNorm <- readRDS("../finalObjects_v3/t_nk_tissue_norm.Rds")
tnkNorm <- tnkNorm[, colnames(tnkNorm) %in% all_meta$cell[all_meta$label == "TNK_tissue"]]

glomNorm <- readRDS("../finalObjects_v3/glom_norm.Rds")
glomNorm <- glomNorm[, colnames(glomNorm) %in% all_meta$cell[all_meta$label == "GLOM"]]

RowMergeSparseMatrices <- function(mat1, mat2) {
  all.mat <- c(list(mat1), mat2)
  all.colnames <- all.rownames <- vector(
    mode = 'list',
    length = length(x = all.mat)
  )
  for (i in seq_along(along.with = all.mat)) {
    if (is.data.frame(x = all.mat[[1]])) {
      all.mat[[i]] <- as.matrix(x = all.mat[[i]])
    }
    all.rownames[[i]] <- rownames(x = all.mat[[i]])
    all.colnames[[i]] <- colnames(x = all.mat[[i]])
  }
  use.cbind <- all(duplicated(x = all.rownames)[2:length(x = all.rownames)])
  if (isTRUE(x = use.cbind)) {
    new.mat <- do.call(what = cbind, args = all.mat)
  } else {
    all.mat <- lapply(X = all.mat, FUN = as, Class = "RsparseMatrix")
    all.names <- unique(x = unlist(x = all.rownames))
    new.mat <- RowMergeMatricesList(
      mat_list = all.mat,
      mat_rownames = all.rownames,
      all_rownames = all.names
    )
    rownames(x = new.mat) <- make.unique(names = all.names)
  }
  colnames(x = new.mat) <- make.unique(names = unlist(x = all.colnames))
  return(new.mat)
}

allNorm <- RowMergeSparseMatrices(myeloidNorm, bNorm)
allNorm <- RowMergeSparseMatrices(allNorm, tnkNorm)
allNorm <- RowMergeSparseMatrices(allNorm, glomNorm)

cellchat <- createCellChat(object = allNorm, meta = all_meta %>% rename(samples = sample), group.by = "label")
rm(allNorm)
gc()

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data

dir.create("cellChat_outputs")
pdf("cellChat_outputs/CellChatDB_categories.pdf", width = 8, height = 6)
  print(showDatabaseCategory(CellChatDB))
dev.off()

CellChatDB.use <- subsetDB(CellChatDB)

# set the used database in the object
cellchat@DB <- CellChatDB.use

ptm = Sys.time()

cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database

cellchat <- identifyOverExpressedGenes(cellchat)

cellchat <- identifyOverExpressedInteractions(cellchat)
#> The number of highly variable ligand-receptor pairs used for signaling inference is 1,965

execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

cellchat <- computeCommunProb(cellchat, type = "triMean", population.size = TRUE)

cellchat <- aggregateNet(cellchat)

saveRDS(cellchat, "cellchat_immuneInteractions.Rds")

execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

ptm = Sys.time()
groupSize <- as.numeric(table(cellchat@idents))

pdf("interactions.pdf", width = 12, height = 6)
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Interaction weights/strength")
dev.off() 

# Create output folder if it doesn't exist
output_dir <- "ligandReceptorOutputs"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Get all pathways in the CellChat object
all_pathways <- cellchat@netP$pathways  # or pathways.show <- names(cellchat@netP$weight) depending on your CellChat version

# Loop over all pathways
for (pw in all_pathways) {
  message("Processing pathway: ", pw)
  
  # 1️⃣ Hierarchy plot
  pdf(file.path(output_dir, paste0(pw, "_hierarchy.pdf")), width = 10, height = 8)
  vertex.receiver <- seq(1,4)  # Adjust as needed
  netVisual_aggregate(cellchat, signaling = pw, vertex.receiver = vertex.receiver)
  dev.off()
  
  # 2️⃣ Circle plot
  pdf(file.path(output_dir, paste0(pw, "_circle.pdf")), width = 10, height = 8)
  par(mfrow = c(1,1))
  netVisual_aggregate(cellchat, signaling = pw, layout = "circle")
  dev.off()
  
  # 3️⃣ Chord diagram (aggregate)
  pdf(file.path(output_dir, paste0(pw, "_chord.pdf")), width = 10, height = 8)
  par(mfrow = c(1,1))
  netVisual_aggregate(cellchat, signaling = pw, layout = "chord")
  dev.off()
  
  # 4️⃣ Heatmap
  pdf(file.path(output_dir, paste0(pw, "_heatmap.pdf")), width = 10, height = 8)
  par(mfrow = c(1,1))
  print(netVisual_heatmap(cellchat, signaling = pw, color.heatmap = "Reds"))
  dev.off()
  
  # 5️⃣ Chord diagram (cell type grouped)
  pdf(file.path(output_dir, paste0(pw, "_chord_byGroup.pdf")), width = 10, height = 8)
  netVisual_chord_cell(cellchat, signaling = pw, 
                       title.name = paste0(pw, " signaling network"))
  dev.off()
}