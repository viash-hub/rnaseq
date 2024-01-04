#!/bin/bash

echo ">>> Testing $meta_functionality_name"

"$meta_executable" \
  --paired false \
  --bam $meta_resources_dir/test.bam \
  --bai $meta_resources_dir/test.bam.bai \
  --get_output_stats true \
  --output_bam test.deduped.bam \
  --output_stats test.umi_dedup.stats

echo ">>> Checking whether output exists"
[ ! -f test.deduped.bam ] && echo "File 'test.deduped.bam' does not exist!" && exit 1
[ ! -f test.umi_dedup.stats ] && echo "File 'test.umi_dedup.stats' does not exist!" && exit 1

echo "All tests succeeded!"
exit 0