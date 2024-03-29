#!/bin/bash

echo ">>> Testing $meta_functionality_name"

"$meta_executable" \
  --genome_fasta "$meta_resources_dir/genome.fasta" \
  --transcriptome_fasta "$meta_resources_dir/transcriptome.fasta" \
  --salmon_index Salmon 

echo ">>> Checking whether output exists"
[ ! -d "Salmon" ] && echo "Salmon index does not exist!" && exit 1
[ -z "$(ls -A 'Salmon')" ] && echo "Salmon index is empty!" && exit 1

echo "All tests succeeded!"
exit 0
