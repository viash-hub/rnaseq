#!/bin/bash

set -eo pipefail 

filename="$(basename -- $par_input)"

junction_annotation.py \
    -i $par_input \
    -r $par_refgene \
    -o "$filename" \
    -m $par_min_intron \
    -q $par_map_qual
2> $par_output_log

mv $filename.junction.bed $par_output_junction_bed
mv $filename.junction.Interact.bed $par_output_junction_interact
mv $filename.junction.xls $par_output_junction_sheet
mv $filename.junction_plot.r $par_output_plot_r
mv $filename.splice_events.pdf $par_output_splice_events_plot
mv $filename.splice_junction.pdf $par_output_splice_junctions_plot
