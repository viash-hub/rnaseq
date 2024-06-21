#!/bin/bash

set -eo pipefail 

prefix=$(openssl rand -hex 8)

read_duplication.py \
    -i $par_input \
    -o $prefix \
    -u $par_read_count_upper_limit \
    -q $par_map_qual 

mv $prefix.DupRate_plot.pdf $par_output_duplication_rate_plot
mv $prefix.DupRate_plot.r $par_output_duplication_rate_plot_r
mv $prefix.pos.DupRate.xls $par_output_duplication_rate_mapping
mv $prefix.seq.DupRate.xls $par_output_duplication_rate_sequence
