#!/bin/bash

echo ">>> Testing $meta_functionality_name"

gunzip "$meta_resources_dir/genes.gtf"
gunzip "$meta_resources_dir/gfp.fa"

"$meta_executable" \
  --fasta "$meta_resources_dir/genome.fasta" \
  --gtf "$meta_resources_dir/genes.gtf" \
  --additional_fasta "$meta_resources_dir/gfp.fa" \
  --biotype gene_biotype \
  --fasta_output genome_gfp.fasta \
  --gtf_output genome_gfp.gtf

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">>> Checking whether output exists"
[ ! -f "genome_gfp.fasta" ] && echo "File 'genome_gfp.fasta' does not exist!" && exit 1
[ ! -s "genome_gfp.fasta" ] && echo "File 'genome_gfp.fasta' is empty!" && exit 1
[ ! -f "genome_gfp.gtf" ] && echo "File 'genome_gfp.gtf' does not exist!" && exit 1
[ ! -s "genome_gfp.gtf" ] && echo "File 'genome_gfp.gtf' is empty!" && exit 1

echo "All tests succeeded!"
exit 0
