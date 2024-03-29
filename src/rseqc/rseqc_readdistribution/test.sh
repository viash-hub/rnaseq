#!/bin/bash

gunzip "$meta_resources_dir/hg19_RefSeq.bed.gz"

# define input and output for script
input_bam="$meta_resources_dir/Pairend_StrandSpecific_51mer_Human_hg19.bam"
input_bed="$meta_resources_dir/hg19_RefSeq.bed"
output="read_distribution.txt"

# run executable and test
echo "> Running $meta_functionality_name"

"$meta_executable" \
    --input "$input_bam" \
    --refgene "$input_bed" \
    --output "$output"

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> Asserting output file was created"
[ ! -f "$output" ] && echo "$output was not created" && exit 1
[ ! -f "$output" ] && echo "$output is empty" && exit 1

exit 0