library(Seurat)
library(sctransform)
library(data.table)
library(tidyverse)

# Parse args
args <- commandArgs(TRUE)
tmp <- args[1] 
n <- args[2] 

# Read files
counts <- fread(paste0(tmp, '/counts.csv')) %>% as.data.frame()
row.names(counts) <- counts[,1]
counts <- counts[,2:ncol(counts)] %>% t() # Transpose for compatibility
meta <- read.csv(paste0(tmp, '/meta.csv'), row.names=1)

# Seurat, SCT workflow
seurat <- CreateSeuratObject(counts=counts, meta.data=meta)
seurat <- SCTransform(
    seurat, 
    vars.to.regress=c("cycle_diff", "nUMIs", "mito_perc"),
    variable.features.n=n
)

# Extract residuals of HVGs and their attributes
residuals <- seurat@assays$SCT@scale.data %>% t() %>% as.data.frame()
df_var <- L@feature.attributes %>% 
    arrange(desc(residual_variance)) %>% 
    select(residual_variance)

# Write
write.csv(residuals, paste0(tmp, '/residuals.csv'))
write.csv(df_var, paste0(tmp, '/df_var.csv'))

