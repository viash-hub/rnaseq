#!/bin/bash

echo ">>> Testing $meta_functionality_name"

echo ">>> Testing paired-end read samples with multiple replicates"
"$meta_executable" \
  --read_1 $meta_resources_dir/SRR6357070_1.fastq.gz\;$meta_resources_dir/SRR6357071_1.fastq.gz \
  --read_2 $meta_resources_dir/SRR6357070_2.fastq.gz\;$meta_resources_dir/SRR6357071_2.fastq.gz \
  --fastq_1 read_1.merged.fastq \
  --fastq_2 read_2.merged.fastq

echo ">>> Checking whether output exists"
[ ! -f "read_1.merged.fastq" ] && echo "Merged read 1 file does not exist!" && exit 1
[ ! -s "read_1.merged.fastq" ] && echo "Merged read 1 file is empty!" && exit 1
[ ! -f "read_2.merged.fastq" ] && echo "Merged read 2 file does not exist!" && exit 1
[ ! -s "read_2.merged.fastq" ] && echo "Merged read 2 file is empty!" && exit 1

echo ">>> Check number of reads"
rep1_1=$(zcat $meta_resources_dir/SRR6357070_1.fastq.gz | echo $((`wc -l`/4)))
rep1_2=$(zcat $meta_resources_dir/SRR6357070_2.fastq.gz | echo $((`wc -l`/4)))
rep2_1=$(zcat $meta_resources_dir/SRR6357071_1.fastq.gz | echo $((`wc -l`/4)))
rep2_2=$(zcat $meta_resources_dir/SRR6357071_2.fastq.gz | echo $((`wc -l`/4)))
merged_1=$(cat read_1.merged.fastq | echo $((`wc -l`/4)))
merged_2=$(cat read_2.merged.fastq | echo $((`wc -l`/4))) 
[[ $(( $rep1_1 + $rep2_1 )) !=  $merged_1 ]] || [[ $(( $rep1_2 + $rep2_2 )) !=  $merged_2 ]] && echo "Concatenation unsuccessful!" && exit 1

rm read_1.merged.fastq read_2.merged.fastq

echo ">>> Testing single-end read samples with multiple replicates"
"$meta_executable" \
  --read_1 $meta_resources_dir/SRR6357070_1.fastq.gz\;$meta_resources_dir/SRR6357071_1.fastq.gz \
  --fastq_1 read_1.merged.fastq 

echo ">>> Checking whether output exists"
[ ! -f "read_1.merged.fastq" ] && echo "Merged read 1 file does not exist!" && exit 1
[ ! -s "read_1.merged.fastq" ] && echo "Merged read 1 file is empty!" && exit 1

echo ">>> Check number of reads"
rep1_1=$(zcat $meta_resources_dir/SRR6357070_1.fastq.gz | echo $((`wc -l`/4)))
rep2_1=$(zcat $meta_resources_dir/SRR6357071_1.fastq.gz | echo $((`wc -l`/4)))
merged_1=$(cat read_1.merged.fastq | echo $((`wc -l`/4)))
[ $(( $rep1_1 + $rep2_1 )) !=  $merged_1 ] && echo "Concatenation unsuccessful!" && exit 1
echo "All tests succeeded!"
exit 0