#!/bin/bash

gunzip "$meta_resources_dir/hg19_RefSeq.bed.gz"

# define input and output for script
input_bam="$meta_resources_dir/Pairend_StrandSpecific_51mer_Human_hg19.bam"
input_bed="$meta_resources_dir/hg19_RefSeq.bed"
output="strandedness.txt"

# run executable and tests
echo "> Running $meta_functionality_name."

"$meta_executable" \
    --input "$input_bam" \
    --refgene "$input_bed" \
    --output "$output"

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> Checking whether output can be found and has content"

[ ! -f "$output" ] && echo "$output is missing" && exit 1
[ ! -s "$output" ] && echo "$output is empty" && exit 1

exit 0