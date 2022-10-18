#/usr/bin/sh

# Set folder
begin=$1

# Set initial folder
cd $begin 

rm -r report.txt runs results_and_plots
rm -r data/removed_cells
rm -r data/step_*