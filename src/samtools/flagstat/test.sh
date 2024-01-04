#!/bin/bash

echo ">>> Testing $meta_functionality_name"

"$meta_executable" \
  --bam $meta_resources_dir/mapt.NA12156.altex.bam \
  --bai $meta_resources_dir/mapt.NA12156.altex.bam.bai \
  --output mapt.NA12156.altex.flagstat

echo ">>> Checking whether output exists"
[ ! -f mapt.NA12156.altex.flagstat ] && echo "File 'mapt.NA12156.altex.flagstat' does not exist!" && exit 1

echo "All tests succeeded!"
exit 0