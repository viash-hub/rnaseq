#!/bin/bash

echo ">>> Testing $meta_functionality_name"

echo ">>> Generating BAM index"
"$meta_executable" \
  --input $meta_resources_dir/mapt.NA12156.altex.bam \
  --output mapt.NA12156.altex.sorted.bam.bai

echo ">>> Check whether output exists"
[ ! -f mapt.NA12156.altex.sorted.bam.bai ] && echo "File 'mapt.NA12156.altex.sorted.bam.bai' does not exist!" && exit 1

echo "All tests succeeded!"
exit 0