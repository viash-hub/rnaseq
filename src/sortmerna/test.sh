#!/bin/bash

echo ">>> Testing $meta_functionality_name"

echo ">>> Testing for paired-end reads"

"$meta_executable" \
    --paired true \
    --input $meta_resources_dir/SRR6357070_1.fastq.gz,$meta_resources_dir/SRR6357070_2.fastq.gz \
    --ribo_database_manifest $meta_resources_dir/rrna-db-defaults.txt \
    --sortmerna_log  SRR6357070_sortmerna.log \
    --fastq_1  SRR6357070_read1.fastq.gz \
    --fastq_2  SRR6357070_read2.fastq.gz

echo ">> Checking if the correct files are present"
[[ ! -f SRR6357070_read1.fastq.gz ]] || [[ ! -f SRR6357070_read2.fastq.gz ]] && echo "Output fastq file missing" && exit 1
[ ! -f SRR6357070_sortmerna.log ] && echo "Output log file missing" && exit 1

echo ">>> Testing for single-end reads"
"$meta_executable" \
    --paired false \
    --input $meta_resources_dir/SRR6357070_1.fastq.gz \
    ---sortmerna_log  SRR6357070_sortmerna.log \
    --fastq_1  SRR6357070_read1.fastq.gz 
    
echo ">> Checking if the correct files are present"
[ ! -f SRR6357070_read1.fastq.gz ] && echo "Output fastq file missing" && exit 1
[ ! -f SRR6357070_sortmerna.log ] && echo "Output log file missing" && exit 1

echo ">>> Test finished successfully"
exit 0
