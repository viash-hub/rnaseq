#!/bin/bash

echo ">>> Testing $meta_functionality_name"

# find $meta_resources_dir/rRNA -type f > rrna-db-defaults.txt

echo ">>> Testing for paired-end reads"
"$meta_executable" \
    --paired true \
    --input $meta_resources_dir/SRR6357070_1.fastq.gz,$meta_resources_dir/SRR6357070_2.fastq.gz \
    --ribo_database_manifest "$meta_resources_dir/rRNA/silva-arc-16s-id95.fasta;$meta_resources_dir/rRNA/silva-euk-18s-id95.fasta" \
    --sortmerna_log SRR6357070_sortmerna.log \
    --fastq_1 SRR6357070_read_1.fastq.gz \
    --fastq_2 SRR6357070_read_2.fastq.gz

echo ">> Checking if the correct files are present"
[[ ! -f "SRR6357070_read_1.fastq.gz" ]] || [[ ! -f "SRR6357070_read_2.fastq.gz" ]] && echo "Output fastq file is missing!" && exit 1
[[ ! -s "SRR6357070_read_1.fastq.gz" ]] || [[ ! -s "SRR6357070_read_2.fastq.gz" ]] && echo "Output fastq file is empty!" && exit 1
[ ! -f "SRR6357070_sortmerna.log" ] && echo "Output log file is missing!" && exit 1
[ ! -s "SRR6357070_sortmerna.log" ] && echo "Output log file is empty!" && exit 1

rm SRR6357070_read_1.fastq.gz SRR6357070_read_2.fastq.gz SRR6357070_sortmerna.log
rm -rf kvdb/

echo ">>> Testing for single-end reads"
"$meta_executable" \
    --paired false \
    --input $meta_resources_dir/SRR6357070_1.fastq.gz \
    --ribo_database_manifest "$meta_resources_dir/rRNA/silva-arc-16s-id95.fasta;$meta_resources_dir/rRNA/silva-euk-18s-id95.fasta" \
    --sortmerna_log SRR6357070_sortmerna.log \
    --fastq_1 SRR6357070_read_1.fastq.gz 
    
echo ">> Checking if the correct files are present"
[ ! -f "SRR6357070_read_1.fastq.gz" ] && echo "Output fastq file is missing!" && exit 1
[ ! -s "SRR6357070_read_1.fastq.gz" ] && echo "Output fastq file is empty!" && exit 1
[ ! -f "SRR6357070_sortmerna.log" ] && echo "Output log file is missing!" && exit 1
[ ! -s "SRR6357070_sortmerna.log" ] && echo "Output log file is empty!" && exit 1

echo ">>> Test finished successfully"
exit 0

