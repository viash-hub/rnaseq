#!/bin/bash

set -exo pipefail 

filename="$(basename -- $par_input)"

if $par_paired; then
    inner_distance.py \
        -i $par_input \
        -r $par_refgene \
        -o $filename \
        -k $par_sample_size \
        -l $par_lower_bound_size \
        -u $par_upper_bound_size \
        -s $par_step_size \
        -q $par_map_qual \
    > stdout.txt

    head -n 2 stdout.txt > $par_output_stats
    
    mv $filename.inner_distance.txt $par_output_dist
    mv $filename.inner_distance_plot.pdf $par_output_plot
    mv $filename.inner_distance_plot.r $par_output_plot_r
    mv $filename.inner_distance_freq.txt $par_output_freq

fi
