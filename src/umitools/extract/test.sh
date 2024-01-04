#!/bin/bash

echo ">>> Testing $meta_functionality_name"

echo ">>> Testing for paired-end reads"
"$meta_executable" \
    --paired true \
    --input $meta_resources_dir/SRR6357070_1.fastq.gz,$meta_resources_dir/SRR6357070_2.fastq.gz \
    --bc_pattern  NNNN,NNNN \
    --umitools_extract_method string \
    --umitools_umi_separator  '_' \
    --umitools_grouping_method  directional \
    --umi_discard_read  0 \
    --fastq_1 SRR6357070_1.umi_extract.fastq.gz \
    --fastq_2 SRR6357070_2.umi_extract.fastq.gz 

echo ">> Checking if the correct files are present"
[[ ! -f SRR6357070_1.umi_extract.fastq.gz ]] || [[ ! -f SRR6357070_2.umi_extract.fastq.gz ]] && echo "Reads file missing" && exit 1

echo ">>> Testing for umi_discard_reads for paired-end reads"
"$meta_executable" \
    --paired true \
    --input $meta_resources_dir/SRR6357070_1.fastq.gz,$meta_resources_dir/SRR6357070_2.fastq.gz \
    --bc_pattern  NNNN,NNNN \
    --umitools_extract_method string \
    --umitools_umi_separator  '_' \
    --umitools_grouping_method  directional \
    --umi_discard_read  2 \
    --fastq_1 SRR6357070_1.umi_extract.fastq.gz \

echo ">> Checking if the correct files are present"
[ ! -f SRR6357070_1.umi_extract.fastq.gz ] && echo "Read 1 file missing" && exit 1
[ -f SRR6357070_1.umi_extract.fastq.gz ] && echo "Read 2 not discarded" && exit 1

echo ">>> Testing for single-end reads"
"$meta_executable" \
    --paired false \
    --input $meta_resources_dir/SRR6357070_1.fastq.gz \
    --bc_pattern  NNNN \
    --umitools_extract_method string \
    --umitools_umi_separator  '_' \
    --umitools_grouping_method  directional \
    --umi_discard_read  0 \
    --fastq_1 SRR6357070_1.umi_extract.fastq.gz 

echo ">> Checking if the correct files are present"
[ ! -f SRR6357070_1.umi_extract.fastq.gz ] && echo "Trimmed reads file missing" && exit 1

echo ">>> Test finished successfully"
exit 0
