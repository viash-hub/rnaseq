#!/bin/bash

gunzip "$meta_resources_dir/genes.gtf.gz"

echo ">>>Testing $metat_functionality_name"
"$meta_executable" \
  --fasta "$meta_resources_dir/genome.fasta" \
  --gtf "$meta_resources_dir/genes.gtf" \
  --filtered_gtf filtered_genes.gtf

echo ">>> Check whether output exists"
[ ! -f "filtered_genes.gtf" ] && echo "Filtered GTF file does not exist!" && exit 1
[ ! -s "filtered_genes.gtf" ] && echo "Filtered GTF file is empty!" && exit 1

echo "All tests succeeded!"
exit 0