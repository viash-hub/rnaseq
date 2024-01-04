#!/bin/bash

echo ">>> Testing $meta_functionality_name"

echo ">>> Generating BAM index"
"$meta_executable" \
  --input $meta_resources_dir/mapt.NA12156.altex.bam \
  --bam_csi_index false \
  --output_bai mapt.NA12156.altex.bam.bai

echo ">>> Check whether output exists"
[ ! -f mapt.NA12156.altex.bam.bai ] && echo "File 'mapt.NA12156.altex.bam.bai' does not exist!" && exit 1

echo ">>> Generating CSI index"
"$meta_executable" \
  --input $meta_resources_dir/mapt.NA12156.altex.bam \
  --bam_csi_index true \
  --output_csi mapt.NA12156.altex.bam.csi

echo ">>> Check whether output exists"
[ ! -f mapt.NA12156.altex.bam.csi ] && echo "File 'mapt.NA12156.altex.bam.csi' does not exist!" && exit 1

echo "All tests succeeded!"
exit 0