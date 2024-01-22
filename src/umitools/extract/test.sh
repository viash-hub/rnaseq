#!/bin/bash

echo ">>> Testing $meta_functionality_name"

echo ">>> Testing for paired-end reads"
"$meta_executable" \
    --paired true \
    --input $meta_resources_dir/scrb_seq_fastq.1.gz,$meta_resources_dir/scrb_seq_fastq.2.gz \
    --bc_pattern CCCCCCNNNNNNNNNN,CCCCCCNNNNNNNNNN \
    --umitools_extract_method string \
    --umitools_umi_separator  '_' \
    --umitools_grouping_method  directional \
    --umi_discard_read  0 \
    --fastq_1 scrb_seq_fastq.1.umi_extract.fastq.gz \
    --fastq_2 scrb_seq_fastq.2.umi_extract.fastq.gz 

echo ">> Checking if the correct files are present"
[[ ! -f scrb_seq_fastq.1.umi_extract.fastq.gz ]] || [[ ! -f scrb_seq_fastq.2.umi_extract.fastq.gz ]] && echo "Reads file missing" && exit 1
[ ! -s "scrb_seq_fastq.1.umi_extract.fastq.gz" ] && echo "Read 1 file is empty" && exit 1
[ ! -s "scrb_seq_fastq.2.umi_extract.fastq.gz" ] && echo "Read 2 file is empty" && exit 1

rm scrb_seq_fastq.1.umi_extract.fastq.gz scrb_seq_fastq.2.umi_extract.fastq.gz

echo ">>> Testing for paired-end reads with umi_discard_reads option"
"$meta_executable" \
    --paired true \
    --input $meta_resources_dir/scrb_seq_fastq.1.gz,$meta_resources_dir/scrb_seq_fastq.2.gz \
    --bc_pattern CCCCCCNNNNNNNNNN,CCCCCCNNNNNNNNNN \
    --umitools_extract_method string \
    --umitools_umi_separator '_' \
    --umitools_grouping_method directional \
    --umi_discard_read 2 \
    --fastq_1 scrb_seq_fastq.1.umi_extract.fastq.gz \

echo ">> Checking if the correct files are present"
[ ! -f "scrb_seq_fastq.1.umi_extract.fastq.gz" ] && echo "Read 1 file is missing" && exit 1
[ ! -s "scrb_seq_fastq.1.umi_extract.fastq.gz" ] && echo "Read 1 file is empty" && exit 1
[ -f "scrb_seq_fastq.2.umi_extract.fastq.gz" ] && echo "Read 2 is not discarded" && exit 1

echo ">>> Testing for single-end reads"
"$meta_executable" \
    --paired false \
    --input $meta_resources_dir/slim.fastq.gz \
    --bc_pattern  "^(?P<umi_1>.{3}).{4}(?P<umi_2>.{2})" \
    --umitools_extract_method regex \
    --umitools_umi_separator '_' \
    --umitools_grouping_method directional \
    --umi_discard_read 0 \
    --fastq_1 slim.umi_extract.fastq.gz 

echo ">> Checking if the correct files are present"
[ ! -f "slim.umi_extract.fastq.gz" ] && echo "Trimmed reads file missing" && exit 1
[ ! -s "slim.umi_extract.fastq.gz" ] && echo "Trimmed reads file is empty" && exit 1

echo ">>> Test finished successfully"
exit 0
