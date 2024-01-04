#!/bin/bash

gunzip "$meta_resources_dir/genes.gtf.gz"

echo ">>> Testing $meta_functionality_name"

"$meta_executable" \
  --salmon_quant_results $meta_resources_dir/WT_REP1.salmon_quant,$meta_resources_dir/WT_REP1.salmon_quant \
  --gtf $meta_resources_dir/genes.gtf \
  --gtf_extra_attributes gene_name \
  --gtf_group_features gene_id \
  --tsv "salmon_tx2gene.tsv" 

echo ">>> Checking whether output exists"
[ ! -f salmon_tx2gene.tsv ] && echo "File 'salmon_tx2gene.tsv' does not exist!" && exit 1

echo "All tests succeeded!"
exit 0
