#!/bin/bash

echo ">>> Testing $meta_functionality_name"

"$meta_executable" \
  --bam "$meta_resources_dir/sample.genome.bam" \
  --fasta "$meta_resources_dir/genome.fasta" \
  --extra_picard_args "--REMOVE_DUPLICATES false" \
  --output_bam "sample.MarkDuplicates.genome.bam" \
  --metrics "sample.MarkDuplicates.metrics.txt"

echo ">>> Check whether output exists"
[ ! -f sample.MarkDuplicates.genome.bam ] && echo "MarkDuplicates output BAM file does not exist!" && exit 1
[ ! -f sample.MarkDuplicates.metrics.txt ] && echo "MarkDuplicates output metrics file does not exist!" && exit 1

echo "All tests succeeded!"
exit 0