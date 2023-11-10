#!/bin/bash

# define input and output for script
input_bam="SRR6357070.bam"
input_bed="genome_gfp.bed"
output="read_distribution.txt"

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
    --output "$tmpdir/$output"

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> asserting output file was created"

[[ ! -f "$tmpdir/$output" ]] && echo "$output was not created" && exit 1

exit 0