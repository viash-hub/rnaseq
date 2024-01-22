#!/bin/bash

gunzip "$meta_resources_dir/hg19_RefSeq.bed.gz"

# define input and output for script
input_bam="$meta_resources_dir/Pairend_StrandSpecific_51mer_Human_hg19.bam"
input_bed="$meta_resources_dir/hg19_RefSeq.bed"

output_stats="inner_distance_stats.txt"
output_dist="inner_distance.txt"
output_plot="inner_distance_plot.pdf"
output_plot_r="inner_distance_plot.r"
output_freq="inner_distance_freq.txt"

# Run executable
echo "> Running $meta_functionality_name"

"$meta_executable" \
    --input $input_bam \
    --refgene $input_bed \
    --output_stats $output_stats \
    --output_dist $output_dist \
    --output_plot $output_plot \
    --output_plot_r $output_plot_r \
    --output_freq $output_freq

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> asserting output has been created for paired read input"

[ ! -f "$output_stats" ] && echo "$output_stats was not created" && exit 1
[ ! -s "$output_stats" ] && echo "$output_stats is empty" && exit 1
[ ! -f "$output_dist" ] && echo "$output_dist was not created" && exit 1
[ ! -s "$output_dist" ] && echo "$output_dist is empty" && exit 1
[ ! -f "$output_plot" ] && echo "$output_plot was not created" && exit 1
[ ! -s "$output_plot" ] && echo "$output_plot is empty" && exit 1
[ ! -f "$output_plot_r" ] && echo "$output_plot_r was not created" && exit 1
[ ! -s "$output_plot_r" ] && echo "$output_plot_r is empty" && exit 1
[ ! -f "$output_freq" ] && echo "$output_freq was not created" && exit 1
[ ! -s "$output_freq" ] && echo "$output_freq is empty" && exit 1

exit 0