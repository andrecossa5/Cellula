# Gene Set Analysis 

########################################################################

# Libraries
library(reticulate)
library(tidyverse)
library(SingleCellExperiment)
library(scater)
library(scran)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(openxlsx)
library(reticulate)

# Source my functions
for (x in list.files('/Users/IEO5505/Desktop/single_cell_package/')) {
    source(paste0('/Users/IEO5505/Desktop/single_cell_package/', x))
}

########################################################################

# Set paths and options
args = commandArgs(trailingOnly=TRUE)
path_main <- args[1]
step <- args[2]
do_msigdb_ <- args[3]
do_markers <- args[4]
do_GMs <- args[5]
do_PCs <- args[6] 
do_DE <- args[7]

path_main <- '/Users/IEO5505/Desktop/AML_GFP_retention/'
do_msigdb_ <- 'no'
do_markers <- 'no'
do_GMs <- 'no'
do_PCs <- 'no' 
do_DE <- 'yes'

path_data <- paste0(path_main, 'data/')
path_results <- paste0(path_main, 'results_and_plots/')

########################################################################

# Get msigdb gene sets and save it in a single df for future use
dbases <- msigdb_downloader()

if (do_msigdb_ == 'yes') {
    for (name in names(dbases)) {
        dbases[[name]] <- dbases[[name]] %>% mutate(Collection = name)
    }
    df <- dplyr::bind_rows(dbases) 
    df <- merge(
                    df,
                    bitr(df$entrez_gene, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db") %>% 
                        rename(entrez_gene=ENTREZID, symbol_gene=SYMBOL),
                    by='entrez_gene',
                    all.x = TRUE
                ) %>% 
                dplyr::select(c(gs_name, entrez_gene, symbol_gene, Collection))

    # Save as .csv
    write.csv(df, paste0(path_data, 'msigdb.csv'))
}

########################################################################

# MARKERs
path_markers <- paste0(path_results, 'clusters_modules/markers/', step, '/')

##--------------------------------------------------------------------##

# ORA

# Load markers
if (do_markers == 'yes') {
    markers <- read.xlsx(paste0(path_markers, 'markers.xlsx'))

    # Prepare them for ORA
    markers.ls <- create_named_list(unique(markers$cluster))
    for (x in names(markers.ls)) { markers.ls[[x]] <- markers %>% filter(cluster == x) %>% pull(genes) }

    # Here we go

    # Navigate to ORA folder
    path_ORA <- paste0(path_results, 'pathways/ORA/')
    setwd(path_ORA)

    # Create markers folder and cd there
    dir.create('markers')
    setwd('markers')

    # For each cluster markers...
    for (x in names(markers.ls)) { 
        dir.create(x)
        setwd(x)
        genes <- markers.ls[[x]]
        do_ORA(genes, dbases, 0.1, n_in=100, n_out=50, minGSSize=50, maxGSSize=500)
        setwd('..')
    }

    ##--------------------------------------------------------------------##

    # GSEA

    # Load markers
    markers <- read.xlsx(paste0(path_markers, 'GSEA_markers.xlsx'))

    # Prepare them for GSEA
    markers.ls <- create_named_list(unique(markers$cluster))
    for (x in names(markers.ls)) { 
        g_ <- markers %>% filter(cluster == x) %>% pull(log2FC) 
        names(g_) <- markers %>% filter(cluster == x) %>% pull(genes) 
        markers.ls[[x]] <- g_
    }

    # Here we go

    # Navigate to GSEA folder
    path_GSEA <- paste0(path_results, 'pathways/GSEA/')
    setwd(path_GSEA)

    # Create markers folder and cd there
    dir.create('markers')
    setwd('markers')

    # For each cluster markers...
    for (x in names(markers.ls)) {  
        dir.create(x)
        setwd(x)
        genes <- markers.ls[[x]]
        do_GSEA(genes, dbases, 0.1, n_out=50, minGSSize=50, maxGSSize=500)
        setwd('..')
    }
}

########################################################################

# MODULES and SIGNATUREs
path_GMs <- paste0(path_results, 'clusters_modules/modules/')

##--------------------------------------------------------------------##

# ORA
if (do_GMs == 'yes') {
    # Load data
    hotspot <- py_load_object(paste0(path_GMs, 'GMs_hotspot.txt'))
    custom <- py_load_object(paste0(path_GMs, 'GMs_custom.txt'))
    sets.ls <- c(hotspot, custom)

    # Here we go

    # Navigate to ORA folder
    path_ORA <- paste0(path_results, 'pathways/ORA/')
    setwd(path_ORA)

    # Create markers folder and cd there
    dir.create('GMs')
    setwd('GMs')

    # For each set...
    for (x in names(sets.ls)) { 
        dir.create(x)
        setwd(x)
        genes <- sets.ls[[x]]
        do_ORA(genes, dbases, 0.1, n_in=100, n_out=50, minGSSize=50, maxGSSize=500)
        setwd('..')
    }

}

########################################################################

# Pseudobulk PCs
path_PCs <- paste0(path_results, 'pp/')

##--------------------------------------------------------------------##

# ORA

# Load markers
if (do_PCs == 'yes') {
    pcs <- read.csv(paste0(path_PCs, 'pseudobulk/loadings.csv'), row.names=1)

    # Prepare them for ORA
    loadings.ls <- create_named_list(colnames(pcs))
    for (x in names(loadings.ls)) { 
        loadings.ls[[x]] <- row.names(pcs %>% arrange(desc(.data[[x]])))[1:100] 
    }

    # Here we go

    # Navigate to ORA folder
    path_ORA <- paste0(path_results, 'pathways/ORA/')
    setwd(path_ORA)

    # Create markers folder and cd there
    dir.create('pseudobulk_pcs_loadings')
    setwd('pseudobulk_pcs_loadings')

    # For each PC top loadings set
    for (x in names(loadings.ls)) { 
        dir.create(x)
        setwd(x)
        genes <- loadings.ls[[x]]
        do_ORA(genes, dbases, 0.1, n_in=100, n_out=50, minGSSize=50, maxGSSize=500)
        setwd('..')
    }

    ##--------------------------------------------------------------------##

    # GSEA

    # Re-format for GSEA
    pcs <- read.csv(paste0(path_PCs, 'pseudobulk/loadings.csv'), row.names=1)
    loadings.ls <- create_named_list(colnames(pcs))
    # Here we go
    for (x in names(loadings.ls)) { 
        g <- pcs %>% arrange(desc(.data[[x]])) %>% pull(.data[[x]])
        names(g) <- row.names(pcs %>% arrange(desc(.data[[x]])))
        loadings.ls[[x]] <- g
    }

    # Here we go

    # Navigate to GSEA folder
    path_GSEA <- paste0(path_results, 'pathways/GSEA/')
    setwd(path_GSEA)

    # Create markers folder and cd there
    dir.create('pseudobulk_pcs_loadings')
    setwd('pseudobulk_pcs_loadings')

    # For each PCs loading...
    for (x in names(loadings.ls)) {  
        dir.create(x)
        setwd(x)
        genes <- loadings.ls[[x]]
        do_GSEA(genes, dbases, 0.1, n_out=50, minGSSize=50, maxGSSize=500)
        setwd('..')
    }
}

########################################################################


# DE
path_DE <- paste0(path_results, 'DE/')
setwd(path_DE)

##--------------------------------------------------------------------##

# ORA

# Load markers
if (do_DE == 'yes') {

    DE.ls <- py_load_object("BIG_DEGs.txt")

    # Unpack
    L <- c()
    for (n in names(DE.ls)) {
        l <- DE.ls[[n]]
        for (m in names(l)) {
            L <- c(L, l[m])
        }
    }

    # Navigate to ORA folder
    path_ORA <- paste0(path_results, 'pathways/ORA/')
    setwd(path_ORA)

    # Create markers folder and cd there
    dir.create('DE')
    setwd('DE')

    # For each PC top loadings set
    for (x in names(L)) { 
        dir.create(x)
        setwd(x)
        genes <- L[[x]]
        if (length(genes) > 0) {
            do_ORA(genes, dbases, 0.1, n_in=100, n_out=50, minGSSize=50, maxGSSize=500)
        }
        setwd('..')
    }
}

########################################################################