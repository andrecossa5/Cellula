#!/bin/bash

# Set general paths
path_main=/Users/IEO5505/Desktop/AML_GFP_retention/

cd $path_main/code

######################################################################

# Step 1

# 1_pp
remove=no
integration_check=no

# 2_clusters_modules_step
step=initial
clustering=yes
do_markers=yes
do_GMs=no
do_dimred=yes

######################################################################

# HERE WE GO

# A --> First exploration of new cells

# Go to main/code folder
cd $path_main/code/

echo -e "Begin $step step...\n" >> trace.txt

# 1_pp
echo -e "begin 1_pp" >> trace.txt
python 1_pp.py $path_main $remove $integration_check
echo -e "end 1_pp\n" >> trace.txt

# 2_clusters_modules_step
echo -e "Begin 2_clusters_modules_step" >> trace.txt
python 2_clusters_modules_step.py $path_main $step $clustering $do_markers $do_GMs $do_dimred
echo -e "End 2_clusters_modules_step\n" >> trace.txt

#######################################################################

# New options

# All of the below 
chosen=leiden_0.4

# 3_clusters_assess
partition_to_remove=none

# 6_GSA
do_msigdb_=yes
do_markers=yes 
do_GMs=no 
do_PCs=no 
do_DE=no

######################################################################

# B --> viz + assessment

# Cepo
# echo -e "Begin cepo" >> trace.txt
# Rscript cepo.r $path_main $chosen 
# echo -e "End cepo\n" >> trace.txt

# # 3_clusters_assess: pt1
# python 3_clusters_assess.py $path_main $step $chosen $partition_to_remove

# 4_DA
#echo -e "Begin DA" >> trace.txt
#Rscript 4_DA.r $path_main $step $chosen 
#echo -e "End DA\n" >> trace.txt
 
# Viz
#echo -e "Begin viz" >> trace.txt
## categorical
#python ./viz/categorical.py $path_main $step $chosen 
#
## cells_embeddings
#python ./viz/cells_embeddings.py $path_main $step $chosen $do_GMs
#
## fishplot
#python ./viz/fishplot.py $path_main $step $chosen 
#
## markers
#python ./viz/markers.py $path_main $step $chosen 
#
## violins (only curated)
#python ./viz/violin.py $path_main $step $chosen $do_GMs
#echo -e "End viz\n" >> trace.txt

# GSA
# echo -e "Begin GSA: only markers" >> trace.txt
# Rscript 6_GSA.r $path_main $step $do_msigdb_ $do_markers $do_GMs $do_PCs $do_DE
# echo -e "End GSA: only markers \n" >> trace.txt

######################################################################

# Decision: remove chosen cluster 10

# New options
partition_to_remove=10

######################################################################

# Remove cells 

# 3_clusters_assess: pt2
echo -e "Write out to remove cells" >> trace.txt
python 3_clusters_assess.py $path_main $step $chosen $partition_to_remove
echo -e "End initial step \n" >> trace.txt

######################################################################





