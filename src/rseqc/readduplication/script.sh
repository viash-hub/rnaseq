#!/bin/bash

set -eo pipefail 

filename="$(basename -- $par_input)"

read_duplication.py \
    -i $par_input \
    -o "$filename" \
    -u $par_read_count_upper_limit \
    -q $par_map_qual 

mv $filename.DupRate_plot.pdf $par_output_duplication_rate_plot
mv $filename.DupRate_plot.r $par_output_duplication_rate_plot_r
mv $filename.pos.DupRate.xls $par_output_duplication_rate_mapping
mv $filename.seq.DupRate.xls $par_output_duplication_rate_sequence