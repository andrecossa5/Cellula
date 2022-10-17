#!/usr/bin/env Rscript

# Rscript for Differential Abundance analysis 

########################################################################

# Libraries
library(tidyverse)
library(zellkonverter)
library(SingleCellExperiment)
library(openxlsx)

# Source my functions
for (x in list.files('/Users/IEO5505/Desktop/single_cell_package/')) {
    source(paste0('/Users/IEO5505/Desktop/single_cell_package/', x))
}

########################################################################

# Set paths and options
args = commandArgs(trailingOnly=TRUE)
path_main <- args[1]
step <- args[2]
chosen <- args[3]

# path_main <- '/Users/IEO5505/Desktop/AML_GFP_retention'
# step <- 'initial'
# chosen <- 'leiden_0.4'

path_data <- paste0(path_main, '/data/')
path_results <- paste0(path_main, '/results_and_plots/clusters_modules/')

# Set out paths (create a {step} folder in clusters and markers, to write to, if necessary)
if (step != 'final') {
    path_clusters <- paste0(path_results, 'clusters/', step, '/')
} else {
    path_clusters <- paste0(path_results, 'clusters_final/')
}

# Load data
sce <- readH5AD(paste0(path_data, 'clustered.h5ad'))

# Extract matrix and cluster labels
meta <- as.data.frame(colData(sce))
clusters <- meta[[chosen]]

########################################################################

# DA: EF

# Prepare directories
setwd(path_clusters)
dir.create(paste0(path_clusters, 'DA/'))
setwd(paste0(path_clusters, 'DA/'))

# Compute and write EFs
write.xlsx(calc_EF(clusters, meta[['sample']]), 'EF_leiden_sample.xlsx', rowNames=TRUE, overwrite=TRUE)
write.xlsx(calc_EF(clusters, meta[['GFP_status']]), 'EF_leiden_GFP_status.xlsx', rowNames=TRUE, overwrite=TRUE)

# Compute and write Differential abundance between conditions
contrasts <- list(
                    c('bulk_d5_tr', 'bulk_d5_un'),
                    c('GFP_high_d5_tr', 'GFP_high_d5_un')
                )
 
# Calculate diff abundance and save
for (l in contrasts) {
    diff_ab(l[1], l[2], chosen, paste0(l[1], '_vs_', l[2]))
} 

########################################################################





