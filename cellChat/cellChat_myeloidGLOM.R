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

myeloid$groupLabels <- case_when(
  myeloid$final_annotation %in% c("M5. GPNMBhigh NUPR1high Macrophage", "M7. SPP1high FABP5high Macrophage",
                                  "M9. MERTKhigh FABP5high Macrophage", "M11. GPMNBhigh NUPR1low Macrophage") ~ "SAMs",
  grepl("Macrophage", myeloid$final_annotation) ~ "Macrophage",
  grepl("DC", myeloid$final_annotation) ~ "DCs",
  grepl("MAST", myeloid$final_annotation) ~ "MAST",
  grepl("Proliferating", myeloid$final_annotation) ~ "Proliferating",
  TRUE ~ "Monocyte"
)

glom <- readRDS("../finalObjects_v3/glom_meta.Rds") %>% select(cell, sample, dataset, Type, final_annotation) %>% mutate(label = "GLOM")
glom$groupLabels <- glom$final_annotation

all_meta <- bind_rows(
  myeloid, glom,
) %>% filter(dataset == "scRNAseq")

rm(myeloid, glom)
gc()

# Load and catenate normalized expression data
myeloidNorm <- readRDS("../finalObjects_v3/myeloid_tissue_norm.Rds")
myeloidNorm <- myeloidNorm[, colnames(myeloidNorm) %in% all_meta$cell[all_meta$label == "Myeloid_tissue"]]

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

allNorm <- RowMergeSparseMatrices(myeloidNorm, glomNorm)

cellchat <- createCellChat(object = allNorm, meta = all_meta %>% rename(samples = sample), group.by = "groupLabels")
rm(allNorm, megaMeta)
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
cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

saveRDS(cellchat, "cellchat_myeloidGLOMInteractions.Rds")

ptm = Sys.time()
groupSize <- as.numeric(table(cellchat@idents))

pdf("interactions_myeloidGLOM.pdf", width = 12, height = 6)
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off() 

# Create output folder if it doesn't exist
output_dir <- "ligandReceptorOutputs_myeloidGLOM"
if (!dir.exists(output_dir)) dir.create(output_dir)

mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  pdf(paste0(output_dir, "/", rownames(mat)[i], "_myeloidGLOM_interactions.pdf"), width = 6, height = 6)
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  dev.off()
}

# pdf("SPP1_interaction_contributions_myeloidGLOM.pdf", width = 6, height = 3)
#   netAnalysis_contribution(cellchat, signaling = "SPP1")  
# dev.off()

pairLR.SPP1 <- extractEnrichedLR(cellchat, signaling = "SPP1", geneLR.return = FALSE)

SPP1_pairs <- netAnalysis_contribution(cellchat, signaling = "SPP1", return.data = TRUE)
pdf("SPP1_interaction_contributions_myeloidGLOM.pdf", width = 4, height = 3)
  ggplot(SPP1_pairs$LR.contribution, aes(x = contribution, y = forcats::fct_rev(name))) + 
    geom_bar(stat = "identity") + 
    theme_classic(base_size = 18)
dev.off()


# Create output folder if it doesn't exist
output_dir <- "SPP1_myeloidGLOM"
if (!dir.exists(output_dir)) dir.create(output_dir)

for (i in 1:nrow(pairLR.SPP1)) {
  pdf(paste0(output_dir, "/", pairLR.SPP1[i, "interaction_name"], "_SPP1_contribution_myeloidGLOM.pdf"), width = 9, height = 9)
    netVisual_individual(cellchat, signaling = "SPP1", pairLR.use = pairLR.SPP1[i, ], label.edge = T, layout = "chord")
  dev.off()
}


# pdf("TENASCIN_interaction_contributions_myeloidGLOM.pdf", width = 6, height = 3)
#   netAnalysis_contribution(cellchat, signaling = "TENASCIN")  
# dev.off()

pairLR.TENASCIN <- extractEnrichedLR(cellchat, signaling = "TENASCIN", geneLR.return = FALSE)
TNC_pairs <- netAnalysis_contribution(cellchat, signaling = "TENASCIN", return.data = TRUE) 

pdf("TNC_interaction_contributions_myeloidGLOM.pdf", width = 5, height = 3)
  ggplot(TNC_pairs$LR.contribution, aes(x = contribution, y = forcats::fct_rev(name))) + 
    geom_bar(stat = "identity") + 
    theme_classic(base_size = 18)
dev.off()


# Create output folder if it doesn't exist
output_dir <- "TENASCIN_myeloidGLOM"
if (!dir.exists(output_dir)) dir.create(output_dir)

for (i in 1:nrow(pairLR.TENASCIN)) {
  pdf(paste0(output_dir, "/", pairLR.TENASCIN[i, "interaction_name"], "_TENASCIN_contribution_myeloidGLOM.pdf"), width = 9, height = 9)
    netVisual_individual(cellchat, signaling = "TENASCIN", pairLR.use = pairLR.TENASCIN[i, ], layout = "chord")
  dev.off()
}
