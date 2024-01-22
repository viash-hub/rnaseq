#!/bin/bash

echo ">>> Testing $meta_functionality_name"

"$meta_executable" \
  --transcript_fasta "$meta_resources_dir/transcriptome.fasta" \
  --output "processed_transcriptome.fasta" 

echo ">>> Check whether output exists"
[ ! -f "processed_transcriptome.fasta" ] && echo "Processed FASTA file does not exist!" && exit 1
[ ! -s "processed_transcriptome.fasta" ] && echo "Processed FASTA file is empty!" && exit 1

echo "All tests succeeded!"
exit 0