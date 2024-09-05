#!/bin/bash

set -eo pipefail 

prefix=$(openssl rand -hex 8)

read_duplication.py \
    -i $par_input \
    -o $prefix \
    -u $par_read_count_upper_limit \
    -q $par_map_qual 

[[ ! -z "$prefix.DupRate_plot.pdf" ]] && mv $prefix.DupRate_plot.pdf $par_output_duplication_rate_plot
[[ ! -z "$prefix.DupRate_plot.r" ]] && mv $prefix.DupRate_plot.r $par_output_duplication_rate_plot_r
[[ ! -z "$prefix.pos.DupRate.xls" ]] && mv $prefix.pos.DupRate.xls $par_output_duplication_rate_mapping
[[ ! -z "$prefix.seq.DupRate.xls" ]] && mv $prefix.seq.DupRate.xls $par_output_duplication_rate_sequence
