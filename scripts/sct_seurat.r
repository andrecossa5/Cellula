library(Seurat)
library(sctransform)
library(data.table)
library(tidyverse)

# Parse args
args <- commandArgs(TRUE)
tmp <- args[1] 
n <- args[2] 

# Read files
# tmp <- '/Users/IEO5505/Desktop/example_cellula/data/tmp'
# n <- 5000

counts <- fread(paste0(tmp, '/counts.csv')) %>% as.data.frame()
row.names(counts) <- counts[,1]
counts <- counts[,2:ncol(counts)] %>% t() # Transpose for compatibility
meta <- read.csv(paste0(tmp, '/meta.csv'), row.names=1)

# Seurat, SCT workflow
seurat <- CreateSeuratObject(counts=counts, meta.data=meta)
seurat <- SCTransform(
    seurat, 
    vars.to.regress=c("cycle_diff", "nUMIs", "mito_perc"),
    variable.features.n=n,
    return.only.var.genes=TRUE
)

# Extract residuals of HVGs and theirz attributes
HVGs <- seurat@assays$SCT@scale.data@var.features[1:n]
X <- seurat@assays$SCT@scale.data[HVGs,] %>% t() %>% as.data.frame()
print(dim(X))
X <- X %>% mutate(
                cells=row.names(X), 
                .before=colnames(X)[1]
            )
df_var <- seurat@assays$SCT@SCTModel.list$model1@feature.attributes
df_var$HVG <- row.names(df_var) %in% HVGs

# Write
fwrite(residuals, paste0(tmp, '/residuals.csv'))
write.csv(df_var, paste0(tmp, '/df_var.csv'))

# RDS
# saveRDS(seurat, paste0(tmp, 'seurat.rds'))