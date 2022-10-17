#!/bin/bash

# Set general paths
path_main=/Users/IEO5505/Desktop/AML_GFP_retention/
# path_singularity=/hpcnfs/scratch/PGP/acossa/singularity/single_cell_1.3.sif

######################################################################

# Step 1

# 1_pp
remove=yes
integration_check=no

# 2_clusters_modules_step
step=1
clustering=yes
do_markers=yes
do_GMs=yes
do_dimred=yes

######################################################################

# HERE WE GO

# A --> First exploration of new cells

# Go to main/code folder
cd $path_main/code/

# echo -e "Begin new round: step $step...\n" >> trace_1.txt
# 
# # 1_pp
# echo -e "begin 1_pp" >> trace_1.txt
# python 1_pp.py $path_main $remove $integration_check
# echo -e "end 1_pp\n" >> trace_1.txt
# 
# # 2_clusters_modules_step
# echo -e "Begin 2_clusters_modules_step" >> trace_1.txt
# python 2_clusters_modules_step.py $path_main $step $clustering $do_markers $do_GMs $do_dimred
# echo -e "End 2_clusters_modules_step\n" >> trace_1.txt

######################################################################

# # New options

# All of the below 
chosen=leiden_0.275

# 3_clusters_assess
partition_to_remove=none

# 6_GSA
do_msigdb_=no
do_markers=yes 
do_GMs=yes
do_PCs=yes
do_DE=no

######################################################################

# B --> viz + assessment

# Cepo
# echo -e "Begin cepo" >> trace_1.txt
# Rscript cepo.r $path_main $chosen 
# echo -e "End cepo\n" >> trace_1.txt
# 
# 3_clusters_assess: pt1
# python 3_clusters_assess.py $path_main $step $chosen $partition_to_remove
 
# # 4_DA
# echo -e "Begin DA" >> trace_1.txt
# Rscript 4_DA.r $path_main $step $chosen 
# echo -e "End DA\n" >> trace_1.txt

# Viz
# echo -e "Begin viz" >> trace_1.txt
# # categorical
# python ./viz/categorical.py $path_main $step $chosen 
# 
# # cells_embeddings
# python ./viz/cells_embeddings.py $path_main $step $chosen $do_GMs
# 
# # fishplot
# python ./viz/fishplot.py $path_main $step $chosen 
# 
# # markers
# # python ./viz/markers.py $path_main $step $chosen 
# 
# # violins (all)
# python ./viz/violin.py $path_main $step $chosen $do_GMs
# echo -e "End viz\n" >> trace_1.txt
#  
#  # Pseudobulk PCs
# echo -e "Begin Pseudobulk PCA" >> trace_1.txt
# python pseudobulk_pca.py $path_main 
# echo -e "End Pseudobulk PCA\n" >> trace_1.txt
 
# GSA (GMs and PCA also)
# echo -e "Begin GSA: markers + GMs + PCs loadings" >> trace_1.txt
# Rscript 6_GSA.r $path_main $step $do_msigdb_ $do_markers $do_GMs $do_PCs $do_DE
# echo -e "End GSA: markers + GMs + PCs + loadings\n" >> trace_1.txt

######################################################################

# Decision: new cells to remove?? 
# If not, go on...


######################################################################
 
# Cell cycle
# echo -e "Begin cc analysis" >> trace_1.txt
# python cell_cycle.py $path_main $chosen
# echo -e "End cc analysis\n" >> trace_1.txt

# Heterogeneity
echo -e "Begin PV calculations" >> trace_1.txt
Rscript 7_heterogeneity.r $path_main $chosen
echo -e "End PV calculations\n" >> trace_1.txt

# DE

# GSA on DE

# Path meta analysis

# Label transfer

# CNVs

####################################################################