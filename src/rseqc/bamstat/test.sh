#!/bin/bash


## VIASH START
meta_executable="target/docker/rseqc/rseqc_bamstat/rseqc_bamstat"
meta_resources_dir="testData/test"
meta_temp_dir="output/bamstat"
## VIASH END

sample="SRR6357070"

echo "> Running $meta_functionality_name."

"$meta_executable" \
    --input "$meta_resources_dir/test/$sample.bam" \
    --output "$meta_temp_dir/mapping_quality.txt"

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> Checking whether output can be found and has content"

[[ ! -f "$meta_temp_dir"/mapping_quality.txt ]] && echo "mapping_quality.txt file missing" && exit 1
[[ -z "$meta_temp_dir"/mapping_quality.txt ]] && echo "Bamstat file is empty" && exit 1

exit 0