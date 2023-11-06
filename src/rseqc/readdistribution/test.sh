#!/bin/bash

echo "> Running $meta_functionality_name."

sample="SRR6357070"
genome="genome_gfp"

"$meta_executable" \
    --input "$meta_resources_dir/test/$sample.bam" \
    --refgene "$meta_resources_dir/test/$genome.bed" \

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> asserting output file was created"

[[! -f "$meta_temp_dir/$par_output" ]] && echo "$par_output was not created" && exit 1

exit 0