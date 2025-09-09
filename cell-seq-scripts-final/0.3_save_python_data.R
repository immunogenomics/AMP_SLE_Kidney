#!/usr/bin/env Rscript

source("/data/srlab1/jmears/srlab/scripts/sc_utils.R")


library(RcppCNPy)

#-------------------------------------------------------

sample = commandArgs(trailingOnly=TRUE)

#-------------------------------------------------------


	data <- readRDS(paste("../cell-seq-output/filtered_data/", sample, "_meta_filtered_500nGene_1000nUMI_3pctnontargetMT.rds", sep = ""))

	npySave(paste("/data/srlab2/jmears/jupyter/SLE_Phase2_Kidney_Analysis/cell-seq-output/python_data/", 
		sample, 
		"_data_filtered_500nGene_1000nUMI_3pctnontargetMT.py",
		sep = ""),
		as.matrix(data))
	rm(data)
	gc()
