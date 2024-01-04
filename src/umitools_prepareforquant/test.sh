#!/bin/bash

echo ">>> Testing $meta_functionality_name"

echo ">>> Generating BAM index"
"$meta_executable" \
  --input $meta_resources_dir/test.bam \
  --output test.transcriptome_sorted.bam \
  --log test.log

echo ">>> Check whether output exists"
[ ! -f test.transcriptome_sorted.bam ] && echo "File 'test.transcriptome_sorted.bam' does not exist!" && exit 1
[ ! -f test.log ] && echo "File 'test.log' does not exist!" && exit 1

echo "All tests succeeded!"
exit 0