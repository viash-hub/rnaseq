#!/bin/bash

gunzip "$meta_resources_dir/hg19_RefSeq.bed.gz"

# define input and output for script
input_bam="$meta_resources_dir/Pairend_StrandSpecific_51mer_Human_hg19.bam"
input_bed="$meta_resources_dir/hg19_RefSeq.bed"

output_plot="junction_saturation_plot.pdf"
output_plot_r="junction_saturation_plot.r"

# run executable and test
echo "> Running $meta_functionality_name"
"$meta_executable" \
    --input "$input_bam" \
    --refgene "$input_bed" \
    --output_plot_r "$output_plot_r" \
    --output_plot "$output_plot"

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> asserting  all output files were created"

[ ! -f "$output_plot_r" ] && echo "$output_plot_r was not created" && exit 1
[ ! -s "$output_plot_r" ] && echo "$output_plot_r is empty" && exit 1
[ ! -f "$output_plot" ] && echo "$output_plot was not created" && exit 1
[ ! -s "$output_plot" ] && echo "$output_plot is empty" && exit 1

exit 0