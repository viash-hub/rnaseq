#!/bin/bash

echo ">>> Testing $meta_functionality_name"

echo ">>> Testing for paired-end reads"
"$meta_executable" \
    --input "$meta_resources_dir/SRR6357070_1.fastq.gz;$meta_resources_dir/SRR6357070_2.fastq.gz" \
    --extra_args  '--record-count 1000000 --seed 1' \
    --output_1  SRR6357070_1.subsampled.fastq.gz \
    --output_2  SRR6357070_2.subsampled.fastq.gz 

echo ">> Checking if the correct files are present"
[ ! -f "SRR6357070_1.subsampled.fastq.gz" ] && echo "Subsampled FASTQ file for read 1 is missing!" && exit 1
[ ! -s "SRR6357070_1.subsampled.fastq.gz" ] && echo "Subsampled FASTQ file is empty!" && exit 1
[ ! -f "SRR6357070_2.subsampled.fastq.gz" ] && echo "Subsampled FASTQ file for read 2 is missing" && exit 1
[ ! -s "SRR6357070_2.subsampled.fastq.gz" ] && echo "Subsampled FASTQ file is empty" && exit 1

rm SRR6357070_1.subsampled.fastq.gz SRR6357070_2.subsampled.fastq.gz

echo ">>> Testing for single-end reads"
"$meta_executable" \
    --input $meta_resources_dir/SRR6357070_1.fastq.gz \
    --extra_args  '--record-count 1000000 --seed 1' \
    --output_1  SRR6357070_1.subsampled.fastq.gz 
    
echo ">> Checking if the correct files are present"
[ ! -f "SRR6357070_1.subsampled.fastq.gz" ] && echo "Subsampled FASTQ file is missing" && exit 1
[ ! -s "SRR6357070_1.subsampled.fastq.gz" ] && echo "Subsampled FASTQ file is empty" && exit 1

echo ">>> Tests finished successfully"
exit 0

