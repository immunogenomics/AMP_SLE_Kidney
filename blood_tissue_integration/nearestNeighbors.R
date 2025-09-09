library(RANN)
library(dplyr)
library(tidyverse)
library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)
name <- as.character(args[1])
theta <- as.numeric(args[2])

numNeigh = 50

# Read in PCs
harmony <- readRDS(paste0(name, "/harmony_", theta, "/", name, "_combined_hPCs.Rds"))
meta <- read.csv(paste0("PCs/", name, "_combinedMeta.csv"), header = TRUE)

# Find nearest neighbors of all
neighbors <- nn2(harmony, treetype = "bd", k = numNeigh)

saveRDS(neighbors, paste0(name, "/harmony_", theta, "/nn.Rds"))

calledCellType <- apply(neighbors$nn.idx, 1:2, function(x) meta$annotation[x])
calledCellType <- cbind(meta$annotation, calledCellType)

tissueOfOrigin <- apply(neighbors$nn.idx, 1:2, function(x, row) meta$origin[x])
tissueOfOrigin <- cbind(meta$annotation, tissueOfOrigin)

counters <- apply(tissueOfOrigin, 1, function(x){
  type <- x[1]
  PBMC <- sum(x[2:ncol(tissueOfOrigin)] == "PBMC")
  Tissue <- (numNeigh - PBMC)
  
  c(type, PBMC, Tissue)
})

counters <- as.data.frame(t(counters))
colnames(counters) <- c("annotation", "PBMC", "tissue")

write.csv(counters, paste0(name, "/harmony_", theta, "/pbmcVsTissue.csv"))

collapsed <- counters %>%
  group_by(annotation) %>%
  summarise(
    PBMC = sum(as.integer(PBMC)),
    tissue = sum(as.integer(tissue))
  )

# Tissue
pbmcs <- as.data.frame(collapsed[collapsed$annotation %in% meta$annotation[meta$origin == "PBMC"], ])
tissue <- as.data.frame(collapsed[collapsed$annotation %in% meta$annotation[meta$origin == "Tissue"], ])

pbmcs <- pbmcs %>%
  mutate(
    total = PBMC + tissue,
    percent_PBMC = (PBMC / total) * 100,
    percent_Tissue = (tissue / total) * 100
  ) %>%
  select(annotation, percent_PBMC, percent_Tissue) %>%
  pivot_longer(cols = c(percent_PBMC, percent_Tissue), names_to = "Type", values_to = "Percent")

tissue <- tissue %>%
  mutate(
    total = PBMC + tissue,
    percent_PBMC = (PBMC / total) * 100,
    percent_Tissue = (tissue / total) * 100
  ) %>%
  select(annotation, percent_PBMC, percent_Tissue) %>%
  pivot_longer(cols = c(percent_PBMC, percent_Tissue), names_to = "Type", values_to = "Percent")

# Create the stacked bar chart
pdf(paste0(name, "/harmony_", theta, "/tissueBars.pdf"), width = 10, height = 10)
    ggplot(tissue, aes(x = annotation, y = Percent, fill = Type)) +
    geom_bar(stat = "identity") +
    labs(title = "Percentage of PBMC and Tissue per Cluster (Tissue)",
        x = "Cluster",
        y = "Percentage") +
    scale_fill_manual(values = c("percent_PBMC" = "blue", "percent_Tissue" = "red"),
                        labels = c("PBMC", "Tissue")) +
    theme_minimal() + 
    theme(
        legend.title = element_text(hjust = -0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

# Create the stacked bar chart
pdf(paste0(name, "/harmony_", theta, "/pbmcBars.pdf"), width = 13, height = 10)
    ggplot(pbmcs, aes(x = annotation, y = Percent, fill = Type)) +
    geom_bar(stat = "identity") +
    labs(title = "Percentage of PBMC and Tissue per Cluster (PBMC)",
        x = "Cluster",
        y = "Percentage") +
    scale_fill_manual(values = c("percent_PBMC" = "blue", "percent_Tissue" = "red"),
                        labels = c("PBMC", "Tissue")) +
    theme_minimal() + 
    theme(
        legend.title = element_text(hjust = -0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

########
selfNonSelf <- apply(calledCellType, 1, function(x){
  type <- x[1]
  self <- sum(x[2:ncol(tissueOfOrigin)] == type)
  nonSelf <- (numNeigh - self)
  
 c(type, self, nonSelf)
})

selfNonSelf <- as.data.frame(t(selfNonSelf))
colnames(selfNonSelf) <- c("annotation", "self", "nonSelf")

write.csv(selfNonSelf, paste0(name, "/harmony_", theta, "/selfNonSelf.csv"))

selfNonSelf_collapsed <- selfNonSelf %>%
  group_by(annotation) %>%
  summarise(
    self = sum(as.integer(self)),
    nonSelf = sum(as.integer(nonSelf))
  )

# Tissue
pbmcs <- as.data.frame(selfNonSelf_collapsed[selfNonSelf_collapsed$annotation %in% meta$annotation[meta$origin == "PBMC"], ])
tissue <- as.data.frame(selfNonSelf_collapsed[selfNonSelf_collapsed$annotation %in% meta$annotation[meta$origin == "Tissue"], ])

pbmcs <- pbmcs %>%
  mutate(
    total = self + nonSelf,
    percent_self = (self / total) * 100,
    percent_nonSelf = (nonSelf / total) * 100
  ) %>%
  select(annotation, percent_self, percent_nonSelf) %>%
  pivot_longer(cols = c(percent_self, percent_nonSelf), names_to = "Type", values_to = "Percent")

tissue <- tissue %>%
  mutate(
    total = self + nonSelf,
    percent_self = (self / total) * 100,
    percent_nonSelf = (nonSelf / total) * 100
  ) %>%
  select(annotation, percent_self, percent_nonSelf) %>%
  pivot_longer(cols = c(percent_self, percent_nonSelf), names_to = "Type", values_to = "Percent")

# Create the stacked bar chart
pdf(paste0(name, "/harmony_", theta, "/pbmcBars_selfNonSelf.pdf"), width = 13, height = 10)
    ggplot(pbmcs, aes(x = annotation, y = Percent, fill = Type)) +
    geom_bar(stat = "identity") +
    labs(title = "Percentage of self and nonSelf neighbors per Cluster (PBMCs)",
        x = "Cluster",
        y = "Percentage") +
    scale_fill_manual(values = c("percent_nonSelf" = "red", "percent_self" = "blue"),
                        labels = c("nonSelf", "self")) +
    theme_minimal() + 
    theme(
        legend.title = element_text(hjust = -0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

# Create the stacked bar chart
pdf(paste0(name, "/harmony_", theta, "/tissueBars_selfNonSelf.pdf"), width = 10, height = 10)
    ggplot(tissue, aes(x = annotation, y = Percent, fill = Type)) +
    geom_bar(stat = "identity") +
    labs(title = "Percentage of self and nonSelf per Cluster (Tissue)",
        x = "Cluster",
        y = "Percentage") +
    scale_fill_manual(values = c("percent_nonSelf" = "blue", "percent_self" = "red"),
                        labels = c("nonSelf", "self")) +
    theme_minimal() + 
    theme(
        legend.title = element_text(hjust = -0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()