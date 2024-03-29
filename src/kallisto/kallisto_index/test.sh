#!/bin/bash

echo ">>> Testing $meta_functionality_name"

"$meta_executable" \
  --transcriptome_fasta "$meta_resources_dir/transcriptome.fasta" \
  --kallisto_index Kallisto 

echo ">>> Checking whether output exists"
[ ! -f "Kallisto" ] && echo "Kallisto index does not exist!" && exit 1
[ ! -s "Kallisto" ] && echo "Kallisto index is empty!" && exit 1

echo "All tests succeeded!"
exit 0
