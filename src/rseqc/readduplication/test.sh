#!/bin/bash

# define input and output for script
input_bam="SRR6357070.bam"

output_duplication_rate_plot_r="duplication_rate_plot.r"
output_duplication_rate_plot="duplication_rate_plot.pdf"

# create temporary directory
tmpdir=$(mktemp -d "$meta_temp_dir/$meta_functionality_name-XXXXXXXX")
function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT

# run executable and test
"$meta_executable" \
    --input "$meta_resources_dir/$input_bam" \
    --output_duplication_rate_plot_r "$tmpdir/$output_duplication_rate_plot_r" \
    --output_duplication_rate_plot "$tmpdir/$output_duplication_rate_plot"

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> asserting  all output files were created"

[[ ! -f "$tmpdir/$output_duplication_rate_plot_r" ]] && echo "$output_duplication_rate_plot_r was not created" && exit 1
[[ ! -f "$tmpdir/$output_duplication_rate_plot" ]] && echo "$output_duplication_rate_plot was not created" && exit 1

exit 0