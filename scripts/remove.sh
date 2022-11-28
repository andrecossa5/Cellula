#/usr/bin/sh

# Set folder
begin=$1

# Set initial folder
cd $begin 

rm -r runs results_and_plots
cd data
ls $apth_main | grep -v "curated_signatures" | xargs rm -r