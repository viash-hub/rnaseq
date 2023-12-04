#!/bin/bash

gunzip "$meta_resources_dir/genes.gff.gz"

echo ">>> Testing $meta_functioniality_name"
"$meta_executable" \
  --input "$meta_resources_dir/genes.gff" \
  --output genes.gtf

echo ">>> Check whether output exists"
[ ! -f genes.gtf ] && echo "GTF file does not exist!" && exit 1

echo "All tests succeeded!"
exit 0