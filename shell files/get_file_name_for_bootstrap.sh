#!/bin/bash

# extract pattern from file name.

file_path=/Users/apple/Desktop/EWCE/Results_bootstrap/Tables/p30/

for filename in $file_path/*.csv; do
#    echo $filename | sed 's/.*p.*0_\(.*\).csv/\1/'
    pattern=$(echo $filename | sed 's/.*p.*0_\(.*\).csv/\1/')
    rscript 250genes_bootstrap_plot.R $pattern
done

echo plot_done!
