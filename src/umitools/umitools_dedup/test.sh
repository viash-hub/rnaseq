#!/bin/bash

echo ">>> Testing $meta_functionality_name"

"$meta_executable" \
  --paired false \
  --bam $meta_resources_dir/chr19.bam \
  --bai $meta_resources_dir/chr19.bam.bai \
  --get_output_stats true \
  --output_bam chr19.deduped.bam \
  --output_stats chr19.umi_dedup.stats

echo ">>> Checking whether output exists"
[ ! -f "chr19.deduped.bam" ] && echo "File 'chr19.deduped.bam' does not exist!" && exit 1
[ ! -s "chr19.deduped.bam" ] && echo "File 'chr19.deduped.bam' is empty!" && exit 1
[ ! -d "chr19.umi_dedup.stats" ] && echo "Directory 'chr19.umi_dedup.stats' does not exist!" && exit 1
[ -z "$(ls -A 'chr19.umi_dedup.stats')" ] && echo "Directory 'chr19.umi_dedup.stats' is empty!" && exit 1

echo "All tests succeeded!"
exit 0