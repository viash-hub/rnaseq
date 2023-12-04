#!/bin/bash

echo ">>> Test $meta_functionality_name"

# cat > bbsplit_fasta_list.txt << HERE
# sarscov2,$meta_resources_dir/sarscov2.fa
# human,$meta_resources_dir/human.fa
# HERE

echo ">>> Building BBSplit index"
"$meta_executable" \
  --id test \
  --primary_ref $meta_resources_dir/genome.fasta \
  --bbsplit_fasta_list bbsplit_fasta_list.txt \
  --only_build_index false \
  --bbsplit_index BBSplit_index 

echo ">>> Check whether output exists"
[ ! -d BBSplit_index ] && echo "BBSplit index does not exist!" && exit 1

echo ">>> Filtering ribosomal RNA reads"

echo ">>> Testing with single-end reads and BBSplit index"
"$meta_executable" \
  --paired false \
  --input "$meta_resources_dir/SRR6357070_1.fastq.gz" \
  --only_build_index false \
  --bbsplit_index BBSplit_index \
  --fastq_1 filtered_SRR6357070_1.fastq.gz

echo ">>> Check whether output exists"
[ ! -f filtered_SRR6357070_1.fastq.gz ] && echo "Filtered reads file does not exist!" && exit 1

echo ">>> Testing with paired-end reads and BBSplit index"
"$meta_executable" \
  --paired true \
  --input "$meta_resources_dir/SRR6357070_1.fastq.gz,$meta_resources_dir/SRR6357070_2.fastq.gz" \
  --only_build_index false \
  --built_bbsplit_index BBSplit_index \
  --fastq_1 filtered_SRR6357070_1.fastq.gz \
  --fastq_2 filtered_SRR6357070_2.fastq.gz

echo ">>> Check whether output exists"
[ ! -f filtered_SRR6357070_1.fastq.gz ] && echo "Filtered read 1 file does not exist!" && exit 1
[ ! -f filtered_SRR6357070_2.fastq.gz ] && echo "Filtered read 2 file does not exist!" && exit 1

echo ">>> Testing with single-end reads and primary/non-primary FASTA files"
"$meta_executable" \
  --paired false \
  --input "$meta_resources_dir/SRR6357070_1.fastq.gz" \
  --only_build_index false \
  --primary_ref "$meta_resources_dir/genome.fasta" \
  --bbsplit_fasta_list bbsplit_fasta_list.txt \
  --fastq_1 filtered_SRR6357070_1.fastq.gz

echo ">>> Check whether output exists"
[ ! -f filtered_SRR6357070_1.fastq.gz ] && echo "Filtered reads file does not exist!" && exit 1

echo ">>> Testing with paired-end reads and primary/non-primary FASTA files"
"$meta_executable" \
  --paired true \
  --input "$meta_resources_dir/SRR6357070_1.fastq.gz,$meta_resources_dir/SRR6357070_2.fastq.gz" \
  --only_build_index false \
  --primary_ref "$meta_resources_dir/genome.fasta" \
  --bbsplit_fasta_list bbsplit_fasta_list.txt \
  --fastq_1 filtered_SRR6357070_1.fastq.gz \
  --fastq_2 filtered_SRR6357070_2.fastq.gz

echo ">>> Check whether output exists"
[ ! -f filtered_SRR6357070_1.fastq.gz ] && echo "Filtered read 1 file does not exist!" && exit 1
[ ! -f filtered_SRR6357070_2.fastq.gz ] && echo "Filtered read 2 file does not exist!" && exit 1

echo "All tests succeeded!"
exit 0