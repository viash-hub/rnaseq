#!/bin/bash

echo ">>> Testing $meta_functionality_name"

echo ">>> Testing paired-end read samples with multiple replicates"
"$meta_executable" \
  --read_1 $meta_resources_dir/SRR6357070_1.fastq.gz\;$meta_resources_dir/SRR6357071_1.fastq.gz \
  --read_2 $meta_resources_dir/SRR6357070_2.fastq.gz\;$meta_resources_dir/SRR6357071_2.fastq.gz \
  --fastq_1 read_1.merged.fastq.gz \
  --fastq_2 read_2.merged.fastq.gz

echo ">>> Checking whether output exists"
[ ! -f "read_1.merged.fastq.gz" ] && echo "Merged read 1 file does not exist!" && exit 1
[ ! -s "read_1.merged.fastq.gz" ] && echo "Merged read 1 file is empty!" && exit 1
[ ! -f "read_2.merged.fastq.gz" ] && echo "Merged read 2 file does not exist!" && exit 1
[ ! -s "read_2.merged.fastq.gz" ] && echo "Merged read 2 file is empty!" && exit 1

rm read_1.merged.fastq.gz read_2.merged.fastq.gz

echo ">>> Testing single-end read samples with multiple replicates"
"$meta_executable" \
  --read_1 $meta_resources_dir/SRR6357070_1.fastq.gz\;$meta_resources_dir/SRR6357071_1.fastq.gz \
  --fastq_1 read_1.merged.fastq.gz 

echo ">>> Checking whether output exists"
[ ! -f "read_1.merged.fastq.gz" ] && echo "Merged read 1 file does not exist!" && exit 1
[ ! -s "read_1.merged.fastq.gz" ] && echo "Merged read 1 file is empty!" && exit 1

echo "All tests succeeded!"
exit 0