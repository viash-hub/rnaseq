#!/bin/bash

echo ">>> Testing $meta_functionality_name"

echo ">>> Testing with BAM input"
"$meta_executable" \
  --paired false \
  --input "$meta_resources_dir/SRR1106616_5M_subset.bam" \
  --output lc_extrap.txt 

echo ">>> Check whether output exists"
[ ! -f "lc_extrap.txt" ] && echo "Output file does not exist!" && exit 1
[ ! -s "lc_extrap.txt" ] && echo "Output file is empty!" && exit 1

rm lc_extrap.txt

echo ">>> Testing with BED input"
"$meta_executable" \
  --paired false \
  --input "$meta_resources_dir/a.sorted.bed" \
  --output lc_extrap.txt 

echo ">>> Check whether output exists"
[ ! -f "lc_extrap.txt" ] && echo "Output file does not exist!" && exit 1
[ ! -s "lc_extrap.txt" ] && echo "Output file is empty!" && exit 1

echo "All tests succeeded!"
exit 0