#!/bin/bash

echo ">>> Testing $meta_functionality_name"

"$meta_executable" \
  --paired true \
  --bam "$meta_resources_dir/genome.bam" \
  --output lc_extrap.txt 

echo ">>> Check whether output exists"
[ ! -f lc_extrap.txt ] && echo "Output file does not exist!" && exit 1

echo "All tests succeeded!"
exit 0