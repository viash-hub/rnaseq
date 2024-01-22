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

# Version
text="${meta_functionality_name}:
    rseqc: $(read_duplication.py --version | sed -e 's/read_duplication.py //g')"

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
    mv "$par_versions" "$par_updated_versions"
else
    echo "$text" > "$par_updated_versions"
fi