#!/bin/bash

echo ">>> Testing $meta_functionality_name"

"$meta_executable" \
  --bam $meta_resources_dir/chr19.bam \
  --bai $meta_resources_dir/chr19.bam.bai \
  --output chr19.idxstats

echo ">>> Checking whether output exists"
[ ! -f "chr19.idxstats" ] && echo "File 'chr19.idxstats' does not exist!" && exit 1
[ ! -s "chr19.idxstats" ] && echo "File 'chr19.idxstats' is empty!" && exit 1

echo "All tests succeeded!"
exit 0