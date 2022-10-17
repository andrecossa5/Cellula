#!/usr/bin/env Rscript

# Viz JIs pathways

########################################################################

# Libraries
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# Paths 
# Set paths and options
args = commandArgs(trailingOnly=TRUE)
path_main <- args[1]

#path_main <- '/Users/IEO5505/Desktop/AML_GFP_retention/'
path_results <- paste0(path_main, 'results_and_plots/pathways/')

# Load matrices
original <- read.csv(paste0(path_results, 'Js_original.csv'), row.names=1) %>% as.matrix()
ours <- read.csv(paste0(path_results, 'Js_our.csv'), row.names=1) %>% as.matrix()

summary(original %>% as.numeric())
summary(ours %>% as.numeric())


# Original 
col_fun = colorRamp2(c(0.01, 0.1, 0.5), c("#0066CC", "white", "#CC6600"))
p <- Heatmap(original, clustering_method_rows = "average", show_row_dend = FALSE,
    clustering_method_columns = "average", show_column_dend = FALSE, 
    row_names_side = 'right', row_names_gp= gpar(fontsize = 2), 
    show_column_names = FALSE, border_gp = gpar(col = "black"), 
    row_split = 15, column_split = 15, column_gap = unit(0.5, "mm"), row_gap = unit(0.7, "mm"),
    row_title_gp = gpar(fontsize = 5),
    column_title_gp = gpar(fontsize = 5), col=col_fun)

pdf(paste0(path_results, 'original.pdf'), width=13.5, height=12)
print(p)
dev.off()

# Ours
p <- Heatmap(ours, col=col_fun, clustering_method_rows = "average", show_row_dend = FALSE,
    clustering_method_columns = "average", show_column_dend = FALSE, 
    row_names_side = 'right', row_names_gp= gpar(fontsize = 2), 
    show_column_names = FALSE, border_gp = gpar(col = "black"), 
    row_split = 15, column_split = 15, column_gap = unit(0.5, "mm"), row_gap = unit(0.7, "mm"),
    row_title_gp = gpar(fontsize = 5),
    column_title_gp = gpar(fontsize = 5))

pdf(paste0(path_results, 'ours.pdf'), width=13.5, height=12)
print(p)
dev.off()

########################################################################