#!/bin/bash

# define input and output for script
input_bam="SRR6357070.bam"
input_bed="genome_gfp.bed"

output_plot="junction_saturation_plot.pdf"
output_plot_r="junction_saturation_plot.r"

# create temporary directory
tmpdir=$(mktemp -d "$meta_temp_dir/$meta_functionality_name-XXXXXXXX")
function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT

# run executable and test
echo "> Running $meta_functionality_name, writing to tmpdir $tmpdir."
"$meta_executable" \
    --input "$meta_resources_dir/$input_bam" \
    --refgene "$meta_resources_dir/$input_bed" \
    --output_plot_r "$tmpdir/$output_plot_r" \
    --output_plot "$tmpdir/$output_plot"

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> asserting  all output files were created"

[[ ! -f "$tmpdir/$output_plot_r" ]] && echo "$output_plot_r was not created" && exit 1
[[ ! -f "$tmpdir/$output_plot" ]] && echo "$output_plot was not created" && exit 1

exit 0