#!/usr/bin/env Rscript

source("/data/srlab1/jmears/srlab/scripts/sc_utils.R")
source("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/utils/plotting.R")

#-------------------------------------------------------

filepath1 <- "/data/srlab2/public/srcollab/AMP_Phase_2/SLE/cellranger-5.0.1/GRCh38/"
filepath2 <- "/outs/raw_feature_bc_matrix/"
samples <- read.delim("/data/srlab2/public/srcollab/AMP_Phase_2/SLE/cellranger-5.0.1/GRCh38/list_of_nuc_samples.txt",
                            sep = "\n", header = F)
samples <- as.character(samples[,1])

#-------------------------------------------------------

## FILTERING PARAMETERS
### NGENE : 500
### MT    : NONE

nGenepar <- 500

#-------------------------------------------------------

# COLUMNS TO ADD TO QC DATAFRAME

Num.droplets.500.genes <- c()
Num.droplets.1000.genes <- c()
Num.droplets.500to1000.genes <- c()
Num.droplets.100to200.genes <- c()
Num.droplets.200to300.genes <- c()
Num.droplets.300to400.genes <- c()
Num.droplets.400to500.genes <- c()
Num.droplets.500to600.genes <- c()
Num.droplets.600to700.genes <- c()
Num.droplets.700to800.genes <- c()
Num.droplets.800to900.genes <- c()
Num.droplets.900to1000.genes <- c()


for (i in 1:length(samples)) {

	raw_meta <- readRDS(paste("../nuc-seq-output/raw_metadata/",
                            samples[i],
                            "_raw_meta.rds",
                            sep = "")
			)

## Read in metadata
    meta <- readRDS(paste("../nuc-seq-output/full_metadata/",
                            samples[i],
                            "_meta_filtered_",
                            nGenepar,
                            "nGene.rds",
                            sep = ""))

# extra columns to add to QC_df
	Num.droplets.500.genes[i] <- nrow(meta)
	Num.droplets.1000.genes[i] <- nrow(meta[meta$nGene > 1000,])
	Num.droplets.500to1000.genes[i] <- nrow(raw_meta[raw_meta$nGene > 499 & raw_meta$nGene < 1001,])
	Num.droplets.100to200.genes[i] <- nrow(raw_meta[raw_meta$nGene > 99 & raw_meta$nGene < 201,])
	Num.droplets.200to300.genes[i] <- nrow(raw_meta[raw_meta$nGene > 199 & raw_meta$nGene < 301,])
	Num.droplets.300to400.genes[i] <- nrow(raw_meta[raw_meta$nGene > 299 & raw_meta$nGene < 401,])
	Num.droplets.400to500.genes[i] <- nrow(raw_meta[raw_meta$nGene > 399 & raw_meta$nGene < 501,])
	Num.droplets.500to600.genes[i] <- nrow(raw_meta[raw_meta$nGene > 499 & raw_meta$nGene < 601,])
	Num.droplets.600to700.genes[i] <- nrow(raw_meta[raw_meta$nGene > 599 & raw_meta$nGene < 701,])
	Num.droplets.700to800.genes[i] <- nrow(raw_meta[raw_meta$nGene > 699 & raw_meta$nGene < 801,])
	Num.droplets.800to900.genes[i] <- nrow(raw_meta[raw_meta$nGene > 799 & raw_meta$nGene < 901,])
	Num.droplets.900to1000.genes[i] <- nrow(raw_meta[raw_meta$nGene > 899 & raw_meta$nGene < 1001,])


# make QC_df

	if (i == 1) {
		QC_df <- read.table(paste("/data/srlab2/public/srcollab/AMP_Phase_2/SLE/cellranger-5.0.1/GRCh38/", samples[i], "/outs/metrics_summary.csv", sep = ""),
                   			sep = ",",
                   			header = T)
	}

	else {

		tmp <- read.table(paste("/data/srlab2/public/srcollab/AMP_Phase_2/SLE/cellranger-5.0.1/GRCh38/", samples[i], "/outs/metrics_summary.csv", sep = ""),
                                        sep = ",",
                                        header = T)
		QC_df <- rbind(QC_df, tmp)
	}




rm(raw_meta)
rm(meta)
gc()

}

QC_df$samples <- samples
QC_df$Num.droplets.500.genes <- Num.droplets.500.genes
QC_df$Num.droplets.1000.genes <- Num.droplets.1000.genes 
QC_df$Num.droplets.500to1000.genes <- Num.droplets.500to1000.genes
QC_df$Num.droplets.100to200.genes <- Num.droplets.100to200.genes
QC_df$Num.droplets.200to300.genes  <- Num.droplets.200to300.genes
QC_df$Num.droplets.300to400.genes <- Num.droplets.300to400.genes
QC_df$Num.droplets.400to500.genes <- Num.droplets.400to500.genes
QC_df$Num.droplets.500to600.genes <- Num.droplets.500to600.genes
QC_df$Num.droplets.600to700.genes <- Num.droplets.600to700.genes
QC_df$Num.droplets.700to800.genes <- Num.droplets.700to800.genes
QC_df$Num.droplets.800to900.genes <- Num.droplets.800to900.genes
QC_df$Num.droplets.900to1000.genes <- Num.droplets.900to1000.genes

saveRDS(QC_df, "/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/nuc-seq-output/QC_output/QC_output_all.rds")

