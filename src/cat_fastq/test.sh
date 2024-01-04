#!/bin/bash

echo ">>> Testing $meta_functionality_name"

echo ">>> Testing paired-end read samples with multiple replicates"
"$meta_executable" \
  --read_1 $meta_resources_dir/read1_replicate1 $meta_resources_dir/read1_replicate2 \
  --read2 $meta_resources_dir/read2_replicate1 $meta_resources_dir/read2_replicate2 \
  --fastq_1 read1.merged.fastq.gz \
  --fastq_2 read2.merged.fastq.gz

echo ">>> Checking whether output exists"
[ ! -f read1.merged.fastq.gz ] && echo "Merged read 1 file does not exist!" && exit 1
[ ! -f read2.merged.fastq.gz ] && echo "Merged read 2 file does not exist!" && exit 1

echo ">>> Testing single-end read samples with multiple replicates"
"$meta_executable" \
  --read_1 $meta_resources_dir/read1_replicate1 $meta_resources_dir/read1_replicate2 \
  --fastq_1 read1.merged.fastq.gz 

echo ">>> Checking whether output exists"
[ ! -f read1.merged.fastq.gz ] && echo "Merged read 1 file does not exist!" && exit 1

echo "All tests succeeded!"
exit 0