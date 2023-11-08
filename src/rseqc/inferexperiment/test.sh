#!/bin/bash

sample="SRR6357070"
genome="genome_gfp"

echo "> Running $meta_functionality_name."

"$meta_executable" \
    --input "$meta_resources_dir/test/$sample.bam" \
    --refgene "$meta_resources_dir/test/$genome.bed" \
    --output "$meta_temp_dir/strandness.txt"

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> Checking whether output can be found and has content"

[[ ! -f "$meta_temp_dir/strandness.txt" ]] && echo "Strandness summary file is missing" && exit 1
[[ -z "$meta_temp_dir/strandness.txt" ]] && echo "Strandness summary file is empty" && exit 1

exit 0