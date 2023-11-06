#!/bin/bash

set -eo pipefail 

filename="$(basename -- $par_input)"

junction_saturation.py \
    -i $par_input \
    -r $par_refgene \
    -o $filename \
    -l $par_sampling_percentile_lower_bound \
    -u $par_sampling_percentile_upper_bound \
    -s $par_sampling_percentile_step \
    -m $par_min_intron \
    -v $par_min_splice_read \
    -q $par_map_qual

mv $filename.junctionSaturation_plot.pdf $par_output_plot
mv $filename.junctionSaturation_plot.r $par_output_plot_r