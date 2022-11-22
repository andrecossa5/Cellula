#/usr/bin/sh

# Set folder
begin=$1

# Set initial folder
cd $begin 

rm -r runs results_and_plots
cd data
rm -v !("curated_signatures")