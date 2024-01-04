#!/bin/bash

echo ">>> Testing $meta_functionality_name"

"$meta_executable" \
  --tpm_gene "$meta_resources_dir/gene_tpm.tsv" \
  --counts_gene "$meta_resources_dir/gene_counts.tsv" \
  --counts_gene_length_scaled "$meta_resources_dir/gene_counts_length_scaled.tsv" \
  --counts_gene_scaled "$meta_resources_dir/gene_counts_scaled.tsv" \
  --tpm_transcript "$meta_resources_dir/transcript_tpm.tsv" \
  --counts_transcript "$meta_resources_dir/transcript_counts.tsv" \
  --tx2gene_tsv "$meta_resources_dir/salmon_tx2gene.tsv" \
  --output "salmon_merged_summarizedexperiment" 

echo ">>> Checking whether output exists"
[ ! -d salmon_merged_summarizedexperiment ] && echo "Directory 'salmon_merged_summarizedexperiment' does not exist!" && exit 1

echo "All tests succeeded!"
exit 0
