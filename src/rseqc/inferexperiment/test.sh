#!/bin/bash

# define input and output for script
input_bam="SRR6357070.bam"
input_bed="genome_gfp.bed"
output="strandedness.txt"

# create temporary directory
tmpdir=$(mktemp -d "$meta_temp_dir/$meta_functionality_name-XXXXXXXX")
function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT

# run executable and tests
echo "> Running $meta_functionality_name."

"$meta_executable" \
    --input "$meta_resources_dir/$input_bam" \
    --refgene "$meta_resources_dir/$input_bed" \
    --output "$tmpdir/$output"

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> Checking whether output can be found and has content"

[[ ! -f "$tmpdir/$output" ]] && echo "$output is missing" && exit 1
[[ -z "$tmpdir/$output" ]] && echo "$output is empty" && exit 1

exit 0