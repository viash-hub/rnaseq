#!/bin/bash

gunzip $meta_resources_dir/genes.gtf.gz
untar -xavf $meta_resources_dir/salmon.tar.gz

echo ">>> Testing $meta_functionality_name for alingment quantification"

"$meta_executable" \
  --paired false \
  --strandedness reverse \
  --input "$meta_resources_dir/SRR6357073_1.fastq.gz" \
  --gtf "$meta_resources_dir/genes.gtf" \
  --transcript_fasta "$meta_resources_dir/transcriptome.fasta" \
  --alignemnt_mode true \
  --extra_salmon_quant_args "" \
  --output "salmon_quant_results" 


echo ">>> Checking whether output exists"
[ ! -d salmon_quant_results ] && echo "Salmon quant output directory does not exist!" && exit 1


echo ">>> Testing $meta_functionality_name for inferring strandedness"

"$meta_executable" \
  --paired false \
  --strandedness auto \
  --input "$meta_resources_dir/SRR6357073_1.fastq.gz" \
  --salmon_index "$meta_resources_dir/salmon" \
  --gtf "$meta_resources_dir/genes.gtf" \
  --transcript_fasta "$meta_resources_dir/transcriptome.fasta" \
  --alignemnt_mode false \
  --lib_type "A" \
  --extra_salmon_quant_args "--skipQuant" \
  --output "salmon_quant_out" \
  --json_info "salmon_meta_info.json" 

echo ">>> Checking whether output exists"
[ ! -d salmon_quant_out ] && echo "Salmon quant output directory does not exist!" && exit 1
[ ! -f salmon_meta_info.json ] && echo "Salmon quant meta information file does not exist!" && exit 1

echo "All tests succeeded!"
exit 0
