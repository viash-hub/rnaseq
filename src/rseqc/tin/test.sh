#!/bin/bash

gunzip "$meta_resources_dir/hg19_RefSeq.bed.gz"

# define input and output for script
input_bam="$meta_resources_dir/Pairend_StrandSpecific_51mer_Human_hg19.bam"
input_bai="$meta_resources_dir/Pairend_StrandSpecific_51mer_Human_hg19.bam.bai"
input_bed="$meta_resources_dir/hg19_RefSeq.bed"

output_tin="tin.xls"
output_tin_summary="tin_summary.txt"

# run executable and test
echo "> Running $meta_functionality_name"

"$meta_executable" \
    --bam_input "$input_bam" \
    --bai_input "$input_bai" \
    --refgene "$input_bed" \
    --output_tin "$output_tin" \
    --output_tin_summary "$output_tin_summary"

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> asserting  all output files were created"

[ ! -f $output_tin ] && echo "$output_tin was not created" && exit 1
[ ! -s $output_tin ] && echo "$output_tin is empty" && exit 1
[ ! -f $output_tin_summary ] && echo "$output_tin_summary was not created" && exit 1
[ ! -s $output_tin_summary ] && echo "$output_tin_summary is empty" && exit 1

exit 0