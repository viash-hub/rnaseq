#!/bin/bash

set -exo pipefail 

prefix=$(openssl rand -hex 8)

inner_distance.py \
    -i $par_input \
    -r $par_refgene \
    -o $prefix \
    -k $par_sample_size \
    -l $par_lower_bound_size \
    -u $par_upper_bound_size \
    -s $par_step_size \
    -q $par_map_qual \
> stdout.txt

head -n 2 stdout.txt > $par_output_stats

[[ -f "$prefix.inner_distance.txt" ]] && mv $prefix.inner_distance.txt $par_output_dist
[[ -f "$prefix.inner_distance_plot.pdf" ]] && mv $prefix.inner_distance_plot.pdf $par_output_plot
[[ -f "$prefix.inner_distance_plot.r" ]] && mv $prefix.inner_distance_plot.r $par_output_plot_r
[[ -f "$prefix.inner_distance_freq.txt" ]] && mv $prefix.inner_distance_freq.txt $par_output_freq
