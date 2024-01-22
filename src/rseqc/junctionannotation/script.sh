#!/bin/bash

set -eo pipefail 

prefix=$(openssl rand -hex 8)

junction_annotation.py \
    -i $par_input \
    -r $par_refgene \
    -o $prefix \
    -m $par_min_intron \
    -q $par_map_qual
> $par_output_log

mv $prefix.junction.bed $par_output_junction_bed
mv $prefix.junction.Interact.bed $par_output_junction_interact
mv $prefix.junction.xls $par_output_junction_sheet
mv $prefix.junction_plot.r $par_output_plot_r
mv $prefix.splice_events.pdf $par_output_splice_events_plot
mv $prefix.splice_junction.pdf $par_output_splice_junctions_plot

# Version
text="${meta_functionality_name}:
    rseqc: $(junction_annotation.py --version | sed -e 's/junction_annotation.py //g')"

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
    mv "$par_versions" "$par_updated_versions"
else
    echo "$text" > "$par_updated_versions"
fi