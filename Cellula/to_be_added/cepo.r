#!/usr/bin/env Rscript

# Rscript for differential stabilty analysis 

########################################################################

# Libraries
library(tidyverse)
library(Cepo)
library(zellkonverter)
library(SingleCellExperiment)
library(openxlsx)
# library(reticulate)
# # use_python("/miniconda3/bin/python3.10")
# use_python("~/miniconda/envs/py_single_cell/bin/python")
# sc <- import("scanpy")

########################################################################

# Set paths and options
args = commandArgs(trailingOnly=TRUE)
path_main <- args[1]
chosen <- args[2]

# path_main <- '/Users/IEO5505/Desktop/AML_GFP_retention/'
# chosen <- 'leiden_0.4'

path_data <- paste0(path_main, '/data/')
path_results <- paste0(path_main, '/results_and_plots/clusters_modules/cepo/')

# Load data
sce <- readH5AD(paste0(path_data, '/normalized_complete.h5ad'))
# adata <- sc$read(paste0(path_data, 'normalized_complete.h5ad'))
# sce <- AnnData2SCE(adata)

# Extract matrix and cluster labels
M <- as.matrix(assay(sce, 'X'))
meta <- as.data.frame(colData(sce))
clusters <- meta[[chosen]]

########################################################################

# Cepo 
ds <- Cepo(M, cellTypes=clusters)
 
# ########################################################################

# Save
df <- as.data.frame(ds$stats) 
colnames(df) <- as.character(levels(meta[[chosen]]))
write.csv(df, paste0(path_results, 'cepo_genes.csv'))

########################################################################