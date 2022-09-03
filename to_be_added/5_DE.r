# edgeR DE

########################################################################

# Libraries
library(reticulate)
library(zellkonverter)
library(tidyverse)
library(SingleCellExperiment)
library(edgeR)
library(openxlsx)

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
sce <- readH5AD(paste0(path_data, 'normalized_complete.h5ad'))
meta <- as.data.frame(colData(sce))
M <- assay(sce, 'X')

########################################################################

# Cluster 2.

# vs itself (d5 vs d0), bulk and GFP

# Bulk
name <- '2_vs_itself_d5_vs_d0_bulk'
print(name)
m_test <- meta %>% filter(leiden_0.275 == '2', day %in% c('0', '5'), GFP_status == 'bulk') %>% 
                mutate(var = case_when(day == '0' ~ 'other', day == '5' ~ 'interesting')) 
m_test$var <- as.factor(m_test$var)
res <- do_DE(M[, row.names(m_test)], m_test, batch_cov=FALSE, var='var', ref='other')
write.xlsx(res[['all_genes']], paste0(name, '.xlsx'), rowNames=TRUE, overwrite=TRUE)

# GFP_high
name <- '2_vs_itself_d5_vs_d0_GFP_high'
print(name)
m_test <- meta %>% filter(leiden_0.275 == '2', day %in% c('0', '5'), GFP_status == 'GFP_high') %>% 
                mutate(var = case_when(day == '0' ~ 'other', day == '5' ~ 'interesting')) 
m_test$var <- as.factor(m_test$var)
res <- do_DE(M[, row.names(m_test)], m_test, FALSE, var='var', ref='other')
write.xlsx(res[['all_genes']], paste0(name, '.xlsx'), rowNames=TRUE, overwrite=TRUE)

##-------------------------------------------------------------------##

# vs 0 1 3 (d5 and d0), bulk and GFP

# Bulk, d0
name <- '2_vs_013_d0_bulk'
print(name)
m_test <- meta %>% filter(leiden_0.275 %in% c('0', '1', '2', '3'), day == '0', GFP_status == 'bulk') %>% 
                mutate(var = case_when(
                    leiden_0.275 %in% c('0', '1', '3') ~ 'other', 
                    leiden_0.275 == '2' ~ 'interesting'
                    )
                ) 
m_test$var <- as.factor(m_test$var)
res <- do_DE(M[, row.names(m_test)], m_test, FALSE, var='var', ref='other')
write.xlsx(res[['all_genes']], paste0(name, '.xlsx'), rowNames=TRUE, overwrite=TRUE)

# Bulk, d5
name <- '2_vs_013_d5_bulk'
print(name)
m_test <- meta %>% filter(leiden_0.275 %in% c('0', '1', '2', '3'), day == '5', GFP_status == 'bulk') %>% 
                mutate(var = case_when(
                    leiden_0.275 %in% c('0', '1', '3') ~ 'other', 
                    leiden_0.275 == '2' ~ 'interesting'
                    )
                ) 
m_test$var <- as.factor(m_test$var)
res <- do_DE(M[, row.names(m_test)], m_test, FALSE, var='var', ref='other')
write.xlsx(res[['all_genes']], paste0(name, '.xlsx'), rowNames=TRUE, overwrite=TRUE)

# GFP_high, d0
name <- '2_vs_013_d0_GFP_high'
print(name)
m_test <- meta %>% filter(leiden_0.275 %in% c('0', '1', '2', '3'), day == '0', GFP_status == 'GFP_high') %>% 
                mutate(var = case_when(
                    leiden_0.275 %in% c('0', '1', '3') ~ 'other', 
                    leiden_0.275 == '2' ~ 'interesting'
                    )
                )
m_test$var <- as.factor(m_test$var) 
res <- do_DE(M[, row.names(m_test)], m_test, FALSE, var='var', ref='other')
write.xlsx(res[['all_genes']], paste0(name, '.xlsx'), rowNames=TRUE, overwrite=TRUE)

# GFP_high, d5
name <- '2_vs_013_d5_GFP_high'
print(name)
m_test <- meta %>% filter(leiden_0.275 %in% c('0', '1', '2', '3'), day == '5', GFP_status == 'GFP_high') %>% 
                mutate(var = case_when(
                    leiden_0.275 %in% c('0', '1', '3') ~ 'other', 
                    leiden_0.275 == '2' ~ 'interesting'
                    )
                ) 
m_test$var <- as.factor(m_test$var)
res <- do_DE(M[, row.names(m_test)], m_test, FALSE, var='var', ref='other')
write.xlsx(res[['all_genes']], paste0(name, '.xlsx'), rowNames=TRUE, overwrite=TRUE)

########################################################################

# Cluster 6 7 8 9.

# vs itself (d5 vs d0), bulk and GFP

# Bulk
name <- '6789_vs_itself_d5_vs_d0_bulk'
print(name)
m_test <- meta %>% filter(leiden_0.275 %in% c('6', '7', '8', '9'), day %in% c('0', '5'), GFP_status == 'bulk') %>% 
                mutate(var = case_when(day == '0' ~ 'other', day == '5' ~ 'interesting'))
m_test$var <- as.factor(m_test$var) 
res <- do_DE(M[, row.names(m_test)], m_test, FALSE, var='var', ref='other')
write.xlsx(res[['all_genes']], paste0(name, '.xlsx'), rowNames=TRUE, overwrite=TRUE)

# GFP_high
name <- '6789_vs_itself_d5_vs_d0_GFP_high'
print(name)
m_test <- meta %>% filter(leiden_0.275 %in% c('6', '7', '8', '9'), day %in% c('0', '5'), GFP_status == 'GFP_high') %>% 
                mutate(var = case_when(day == '0' ~ 'other', day == '5' ~ 'interesting')) 
m_test$var <- as.factor(m_test$var)
res <- do_DE(M[, row.names(m_test)], m_test, FALSE, var='var', ref='other')
write.xlsx(res[['all_genes']], paste0(name, '.xlsx'), rowNames=TRUE, overwrite=TRUE)

##-------------------------------------------------------------------##

# vs 0 1 3 (d5 and d0), bulk and GFP

# Bulk, d0
name <- '6789_vs_013_d0_bulk'
print(name)
m_test <- meta %>% filter(leiden_0.275 %in% c('0', '1', '3', '6', '7', '8', '9'), day == '0', GFP_status == 'bulk') %>% 
                mutate(var = case_when(
                    leiden_0.275 %in% c('0', '1', '3') ~ 'other', 
                    leiden_0.275 %in% c( '6', '7', '8', '9') ~ 'interesting'
                    )
                ) 
m_test$var <- as.factor(m_test$var)
res <- do_DE(M[, row.names(m_test)], m_test, FALSE, var='var', ref='other')
write.xlsx(res[['all_genes']], paste0(name, '.xlsx'), rowNames=TRUE, overwrite=TRUE)

# Bulk, d5
name <- '6789_vs_013_d5_bulk'
print(name)
m_test <- meta %>% filter(leiden_0.275 %in% c('0', '1', '3', '6', '7', '8', '9'), day == '5', GFP_status == 'bulk') %>% 
                mutate(var = case_when(
                    leiden_0.275 %in% c('0', '1', '3') ~ 'other', 
                    leiden_0.275 %in% c( '6', '7', '8', '9') ~ 'interesting'
                    )
                ) 
m_test$var <- as.factor(m_test$var)
res <- do_DE(M[, row.names(m_test)], m_test, FALSE, var='var', ref='other')
write.xlsx(res[['all_genes']], paste0(name, '.xlsx'), rowNames=TRUE, overwrite=TRUE)

# GFP_high, d0
name <- '6789_vs_013_d0_GFP_high'
print(name)
m_test <- meta %>% filter(leiden_0.275 %in% c('0', '1', '3', '6', '7', '8', '9'), day == '0', GFP_status == 'GFP_high') %>% 
                mutate(var = case_when(
                    leiden_0.275 %in% c('0', '1', '3') ~ 'other', 
                    leiden_0.275 %in% c( '6', '7', '8', '9') ~ 'interesting'
                    )
                ) 
m_test$var <- as.factor(m_test$var)
res <- do_DE(M[, row.names(m_test)], m_test, FALSE, var='var', ref='other')
write.xlsx(res[['all_genes']], paste0(name, '.xlsx'), overwrite=TRUE)

# GFP_high, d5
name <- '6789_vs_013_d5_GFP_high'
print(name)
m_test <- meta %>% filter(leiden_0.275 %in% c('0', '1', '3', '6', '7', '8', '9'), day == '5', GFP_status == 'GFP_high') %>% 
                mutate(var = case_when(
                    leiden_0.275 %in% c('0', '1', '3') ~ 'other', 
                    leiden_0.275 %in% c( '6', '7', '8', '9') ~ 'interesting'
                    )
                ) 
m_test$var <- as.factor(m_test$var)
res <- do_DE(M[, row.names(m_test)], m_test, FALSE, var='var', ref='other')
write.xlsx(res[['all_genes']], paste0(name, '.xlsx'), overwrite=TRUE)

########################################################################