#!/bin/bash

gunzip "$meta_resources_dir/genes.gtf.gz"

echo ">>> Testing $meta_functionality_name"

"$meta_executable" \
  --salmon_quant_results $meta_resources_dir/WT_REP1.salmon_quant,$meta_resources_dir/WT_REP1.salmon_quant \
  --tx2gene_tsv "$meta_resources_dir/salmon_tx2gene.tsv" \
  --tpm_gene "gene_tpm.tsv" \
  --counts_gene "gene_counts.tsv" \
  --counts_gene_length_scaled "gene_counts_length_scaled.tsv" \
  --counts_gene_scaled "gene_counts_scaled.tsv" \
  --tpm_transcript "transcript_tpm.tsv" \
  --counts_transcript "transcript_counts.tsv"

echo ">>> Checking whether output exists"
[ ! -f gene_tpm.tsv ] && echo "File 'gene_tpm.tsv' does not exist!" && exit 1
[ ! -f gene_counts.tsv ] && echo "File 'gene_counts.tsv' does not exist!" && exit 1
[ ! -f gene_counts_length_scaled.tsv ] && echo "File 'gene_counts_length_scaled.tsv' does not exist!" && exit 1
[ ! -f gene_counts_scaled.tsv ] && echo "File 'gene_counts_scaled.tsv' does not exist!" && exit 1
[ ! -f transcript_tpm.tsv ] && echo "File 'transcript_tpm.tsv' does not exist!" && exit 1
[ ! -f transcript_counts.tsv ] && echo "File 'transcript_counts.tsv' does not exist!" && exit 1

echo "All tests succeeded!"
exit 0
