#!/usr/bin/env Rscript

# Rscript for heterogeneity analysis

########################################################################

# Libraries
library(tidyverse)
library(zellkonverter)
library(SingleCellExperiment)
library(openxlsx)
library(irlba)
library(Rfast)

# Source my functions
for (x in list.files('/Users/IEO5505/Desktop/AML_GFP_retention/single_cell_package/')) {
    source(paste0('/Users/IEO5505/Desktop/AML_GFP_retention/ingle_cell_package/', x))
}


# Utilities

# PV: calculate the phenotypic volume of a certain group of cells
PV <- function(M, cells) {

    M_c <- M[cells, ]
    M_c <- scale(M_c)
    M_c[is.na(M_c)] <- 0
    cov_c <- cova(M_c)
    SVD <- irlba(cov_c, 50)
    pv <- sum(SVD$d^2) / ncol(M_c)

    return(pv)

}


##


# Subsampling routine
subsample <- function(meta, var, x, n_sample, times) {

    cells <- meta[meta[[var]] == x, ] %>% row.names()
    L <- sapply(1:times, function(i) { sample(cells, n_sample) } )

    return (L)

}


################################################################################

# Set paths and options
args = commandArgs(trailingOnly=TRUE)
path_main <- args[1]
chosen <- args[2]

# path_main <- '/hpcnfs/scratch/PGP/acossa/AML_Sara/'
# chosen <- 'leiden_0.275'

path_data <- paste0(path_main, '/data/')
path_results <- paste0(path_main, '/results_and_plots/heterogeneity/')

# Load data
sce <- readH5AD(paste0(path_data, 'clustered.h5ad'))
M <- t(assay(sce, 'X'))

# Extract matrix and cluster labels
meta <- as.data.frame(colData(sce))
clusters <- meta[[chosen]]

##----------------------------------------------------------------------------##

# PV calculations 

# Clusters
var <- chosen
groups <- levels(meta[[var]])
n_samples <- round(min( meta %>% group_by(.data[[var]]) %>% tally() %>% pull(n) ) * 0.8)

# Calculate bootstrapped PVs for each clonal subset
n_times <- 10

PV_matrix <- matrix(0, length(groups), n_times, 
        dimnames=list(groups, paste0('iteration_', as.character(1:n_times))))

# Here we go
for (x in groups) {
    S <- subsample(meta, var, x, n_samples, n_times)
    PV_matrix[x, ] <- sapply( 1:n_times, function(i) { PV(M, S[, i]) } )
    print(paste0('Finished ', as.character(x)))
}

# Write matrix
write.csv(PV_matrix, paste0(path_results, 'PVs_clusters.csv'), row.names=TRUE)

# Read and check
df <- read.csv(paste0(path_results, 'PVs_clusters.csv'), row.names=1)
df$mean <- rowMeans(df)
df <- df %>% arrange(mean)
print(df)

##----------------------------------------------------------------------------##

# Samples
var <- 'sample'
groups <- levels(meta[[var]])
n_samples <- round(min( meta %>% group_by(.data[[var]]) %>% tally() %>% pull(n) ) * 0.8)

# Calculate bootstrapped PVs for each clonal subset
n_times <- 10

PV_matrix <- matrix(0, length(groups), n_times, 
        dimnames=list(groups, paste0('iteration_', as.character(1:n_times))))

# Here we go
for (x in groups) {
    S <- subsample(meta, var, x, n_samples, n_times)
    PV_matrix[x, ] <- sapply( 1:n_times, function(i) { PV(M, S[, i]) } )
    print(paste0('Finished ', as.character(x)))
}

# Write matrix
write.csv(PV_matrix, paste0(path_results, 'PVs_samples.csv'), row.names=TRUE)

# Read and check
df <- read.csv(paste0(path_results, 'PVs_samples.csv'), row.names=1)
df$mean <- rowMeans(df)
df <- df %>% arrange(mean)
print(df)

################################################################################

