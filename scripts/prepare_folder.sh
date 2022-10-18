#/usr/bin/sh

# Set folder
begin=$1

# Set initial folder
cd $begin 

# Create a report.txt file and update data with the removed_cells folder
touch report.txt
cd data 
mkdir removed_cells 

# Go back to main folder, create the results_and_plots and runs folders
cd ..
mkdir results_and_plots runs 

# Update results and plots
cd results_and_plots
mkdir pp clustering vizualization dist_features signatures

# Update vizualization
cd vizualization
mkdir pp clustering dist_features signatures QC

# # Trajectories
# cd ../trajectories
# mkdir WOT cellrank
# cd WOT
# mkdir tmaps trends trajectories transitions
