#!/bin/bash

echo ">>> Testing $meta_functionality_name"

echo ">>> Generating BAM index"
"$meta_executable" \
  --input $meta_resources_dir/chr19.bam \
  --output chr19.sorted.bam.bai

echo ">>> Check whether output exists"
[ ! -f "chr19.sorted.bam.bai" ] && echo "File 'chr19.sorted.bam.bai' does not exist!" && exit 1
[ ! -s "chr19.sorted.bam.bai" ] && echo "File 'chr19.sorted.bam.bai' is empty!" && exit 1

echo "All tests succeeded!"
exit 0