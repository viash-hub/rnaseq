#!/bin/bash

echo "> Running $meta_functionality_name."

sample="SRR6357070"
genome="genome_gfp"

output_plot="junction_saturation_plot.pdf"
output_plot_r="junction_saturation_plot.r"

"$meta_executable" \
    --input "$meta_resources_dir/test/$sample.bam" \
    --refgene "$meta_resources_dir/test/$genome.bed" \
    --output_plot_r "$meta_temp_dir/$output_plot_r" \
    --output_plot "$meta_temp_dir/$output_plot"

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> asserting  all output files were created"

[[ ! -f "$meta_temp_dir/$output_plot_r" ]] && echo "$output_plot_r was not created" && exit 1
[[ ! -f "$meta_temp_dir/$output_plot" ]] && echo "$output_plot was not created" && exit 1

exit 0