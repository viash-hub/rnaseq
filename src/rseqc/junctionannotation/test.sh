#!/bin/bash

gunzip "$meta_resources_dir/hg19_RefSeq.bed.gz"

# define input and output for script
input_bam="$meta_resources_dir/Pairend_StrandSpecific_51mer_Human_hg19.bam"
input_bed="$meta_resources_dir/hg19_RefSeq.bed"

output_junction_bed="junction_annotation.bed"
output_junction_interact="junction_annotation.Interact.bed"
output_junction_sheet="junction_annotation.xls"
output_plot_r="unction_annotation_plot.r"
output_splice_events_plot="splice_events.pdf"
output_splice_junctions_plot="splice_junctions_plot.pdf"
output_log="junction_annotation.log"

# run executable and test
echo "> Running $meta_functionality_name"

"$meta_executable" \
    --input "$input_bam" \
    --refgene "$input_bed" \
    --output_log "$output_log" \
    --output_plot_r "$output_plot_r" \
    --output_junction_bed "$output_junction_bed" \
    --output_junction_interact "$output_junction_interact" \
    --output_junction_sheet "$output_junction_sheet" \
    --output_splice_events_plot "$output_splice_events_plot" \
    --output_splice_junctions_plot "$output_splice_junctions_plot" 

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> asserting  all output files were created"

[ ! -f "$output_log" ] && echo "$output_log was not created" && exit 1
[ ! -f "$output_junction_bed" ] && echo "$output_junction_bed was not created" && exit 1
[ ! -s "$output_junction_bed" ] && echo "$output_junction_bed is empty" && exit 1
[ ! -f "$output_junction_interact" ] && echo "$output_junction_interact was not created" && exit 1
[ ! -s "$output_junction_interact" ] && echo "$output_junction_interact is empty" && exit 1
[ ! -f "$output_junction_sheet" ] && echo "$output_junction_sheet was not created" && exit 1
[ ! -s "$output_junction_sheet" ] && echo "$output_junction_sheet is empty" && exit 1
[ ! -f "$output_plot_r" ] && echo "$output_plot_r was not created" && exit 1
[ ! -s "$output_plot_r" ] && echo "$output_plot_r is empty" && exit 1
[ ! -f "$output_splice_events_plot" ] && echo "$output_splice_events_plot was not created" && exit 1
[ ! -s "$output_splice_events_plot" ] && echo "$output_splice_events_plot is empty" && exit 1
[ ! -f "$output_splice_junctions_plot" ] && echo "$output_splice_junctions_plot was not created" && exit 1
[ ! -s "$output_splice_junctions_plot" ] && echo "$output_splice_junctions_plot is empty" && exit 1

exit 0