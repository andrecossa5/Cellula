# Copycat CNVs

########################################################################

# Libraries
library(reticulate)
library(zellkonverter)
library(tidyverse)
library(SingleCellExperiment)
library(openxlsx)
library(copykat)

# Source my functions
for (x in list.files('/Users/IEO5505/Desktop/single_cell_package/')) {
    source(paste0('/Users/IEO5505/Desktop/single_cell_package/', x))
}

########################################################################

# Set paths
path_main <- '/Users/IEO5505/Desktop/AML_GFP_retention/'
path_data <- paste0(path_main, 'data/')
path_results <- paste0(path_main, 'results_and_plots/DE/')

# Set wd
setwd(path_results)

# Load data
sce <- readH5AD(paste0(path_data, 'raw_matrix.h5ad'))
meta <- as.data.frame(colData(sce))
M <- assay(sce, 'X')

copykat <- copykat(rawmat=as.matrix(M), 
                id.type="S",
                ngene.chr=5, 
                win.size=25, KS.cut=0.1, 
                distance="euclidean", 
                norm.cell.names="",
                genome="hg20",
                n.cores=6)
