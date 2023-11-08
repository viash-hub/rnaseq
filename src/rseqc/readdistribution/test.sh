#!/bin/bash

echo "> Running $meta_functionality_name."

sample="SRR6357070"
genome="genome_gfp"
output="read_distribution.txt"

"$meta_executable" \
    --input "$meta_resources_dir/test/$sample.bam" \
    --refgene "$meta_resources_dir/test/$genome.bed" \
    --output "$meta_temp_dir/$output"

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> asserting output file was created"

[[ ! -f "$meta_temp_dir/$output" ]] && echo "$output was not created" && exit 1

exit 0