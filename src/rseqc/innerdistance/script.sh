#!/bin/bash

set -exo pipefail 

filename="$(basename -- $par_input)"

if $par_paired; then
    inner_distance.py -i $par_input -r $par_refgene -o $filename > stdout.txt
    head -n 2 stdout.txt > $filename.inner_distance_stats.txt
    mv $filename.inner_distance.txt $par_output_dist
    mv $filename.inner_distance_plot.pdf $par_output_plot
    mv $filename.inner_distance_plot.r $par_output_plot_r
    mv $filename.inner_distance_freq.txt $par_output_freq
fi
