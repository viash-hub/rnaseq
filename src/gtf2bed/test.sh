#!/bin/bash

gunzip "$meta_resources_dir/genes.gtf.gz"

echo ">>> Testing $meta_functionality_name"
"$meta_executable" \
  --gtf "$meta_resources_dir/genes.gtf" \
  --bed_output genes.bed

echo ">>> Check whether output exists"
[ ! -f "genes.bed" ] && echo "BED output file does not exist!" && exit 1
[ ! -s "genes.bed" ] && echo "BED output file is empty!" && exit 1

echo "All tests succeeded!"
exit 0