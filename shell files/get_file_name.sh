#!/bin/bash

# extract pattern from file name.

file_path=/Users/apple/Desktop/EWCE/Results_wt_g93a/Tables/p30/

for filename in $file_path/*.csv; do
#    echo $filename | sed 's/.*p.*0_\(.*\).csv/\1/'
    pattern=$(echo $filename | sed 's/.*p.*0_\(.*\).csv/\1/')
    rscript timeplot_11.R $pattern
done

echo plot_done!
