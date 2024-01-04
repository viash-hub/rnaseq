#!/bin/bash

echo ">>> Testing $meta_functionality_name"

"$meta_executable" \
  --bam $meta_resources_dir/mapt.NA12156.altex.bam \
  --bai $meta_resources_dir/mapt.NA12156.altex.bam.bai \
  --output mapt.NA12156.altex.stats

echo ">>> Checking whether output exists"
[ ! -f mapt.NA12156.altex.stats ] && echo "File 'mapt.NA12156.altex.stats' does not exist!" && exit 1

echo "All tests succeeded!"
exit 0