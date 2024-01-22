#!/bin/bash

echo ">>> Testing $meta_functionality_name"

echo ">>> Generating BAM index"
"$meta_executable" \
  --input $meta_resources_dir/chr19.bam \
  --bam_csi_index false \
  --output_bai chr19.bam.bai

echo ">>> Check whether output exists"
[ ! -f chr19.bam.bai ] && echo "File 'mapt.NA12156.altex.bam.bai' does not exist!" && exit 1
[ ! -s chr19.bam.bai ] && echo "File 'mapt.NA12156.altex.bam.bai' is empty!" && exit 1

echo ">>> Generating CSI index"
"$meta_executable" \
  --input $meta_resources_dir/chr19.bam \
  --bam_csi_index true \
  --output_csi chr19.bam.csi

echo ">>> Check whether output exists"
[ ! -f "chr19.bam.csi" ] && echo "File 'mapt.NA12156.altex.bam.csi' does not exist!" && exit 1
[ ! -s "chr19.bam.csi" ] && echo "File 'mapt.NA12156.altex.bam.csi' is empty!" && exit 1

echo "All tests succeeded!"
exit 0