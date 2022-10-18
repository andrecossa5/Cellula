#/usr/bin/sh

# Default workflow. scanpy, 2000 HVGs, red_s:original, 30PCs, 15 NN, 0.2-0.6 res range.
path_main=/Users/IE05505/Desktop/sc_pipeline_prova/

python 1_pp.py -p $path_main --step 0 --norm scanpy --n_HVGs 2000 --score scanpy
python 3_integration_diagnostics.py -p $path_main --step 0 --chosen red_s:original 
python 4_clustering.py -p $path_main --step 0 --range 0.2:0.6 
python 5_clustering_diagnostics.py -p $path_main --step 0  
python 5_clustering_diagnostics.py -p $path_main --step 0 --chosen ...
python 6_signatures.py -p $path_main --step 0 --Hotspot --barkley --wu --scoring scanpy
python 7_dist_features.py -p $path_main --step 0 
# python 8_label_transfer.py -p $path_main --step 0 

