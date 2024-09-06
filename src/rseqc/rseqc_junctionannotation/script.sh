#!/bin/bash

set -eo pipefail 

prefix=$(openssl rand -hex 8)
input="testData/unit_test_resources/test.paired_end.sorted.bam"
refgene="testData/unit_test_resources/test.bed"
junction_annotation.py \
    -i $par_input \
    -r $par_refgene \
    -o $prefix \
    -m $par_min_intron \
    -q $par_map_qual > $par_output_log

[[ -f "$prefix.junction.bed" ]] && mv $prefix.junction.bed $par_output_junction_bed
[[ -f "$prefix.junction.Interact.bed" ]] && mv $prefix.junction.Interact.bed $par_output_junction_interact
[[ -f "$prefix.junction.xls" ]] && mv $prefix.junction.xls $par_output_junction_sheet
[[ -f "$prefix.junction_plot.r" ]] && mv $prefix.junction_plot.r $par_output_plot_r
[[ -f "$prefix.splice_events.pdf" ]] && mv $prefix.splice_events.pdf $par_output_splice_events_plot
[[ -f "$prefix.splice_junction.pdf" ]] && mv $prefix.splice_junction.pdf $par_output_splice_junctions_plot
