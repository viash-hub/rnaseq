#!/bin/bash

# define input and output for script
input_bam="$meta_resources_dir/test.paired_end.sorted.bam"

output_duplication_rate_plot_r="duplication_rate_plot.r"
output_duplication_rate_plot="duplication_rate_plot.pdf"

# run executable and test
echo "> Running $meta_functionality_name"
"$meta_executable" \
    --input "$input_bam" \
    --output_duplication_rate_plot_r "$output_duplication_rate_plot_r" \
    --output_duplication_rate_plot "$output_duplication_rate_plot"

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> asserting  all output files were created"
[ ! -f "$output_duplication_rate_plot_r" ] && echo "$output_duplication_rate_plot_r was not created" && exit 1
[ ! -s "$output_duplication_rate_plot_r" ] && echo "$output_duplication_rate_plot_r is empty" && exit 1
[ ! -f "$output_duplication_rate_plot" ] && echo "$output_duplication_rate_plot was not created" && exit 1
[ ! -s "$output_duplication_rate_plot" ] && echo "$output_duplication_rate_plot is empty" && exit 1

exit 0