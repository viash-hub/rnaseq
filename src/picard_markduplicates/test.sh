#!/bin/bash

echo ">>> Testing $meta_functionality_name"

"$meta_executable" \
  --bam "$meta_resources_dir/a.bam" \
  --fasta "$meta_resources_dir/genome.fasta" \
  --extra_picard_args "--REMOVE_DUPLICATES false" \
  --output_bam "a.MarkDuplicates.genome.bam" \
  --metrics "a.MarkDuplicates.metrics.txt"

echo ">>> Check whether output exists"
[ ! -f "a.MarkDuplicates.genome.bam" ] && echo "MarkDuplicates output BAM file does not exist!" && exit 1
[ ! -s "a.MarkDuplicates.genome.bam" ] && echo "MarkDuplicates output BAM file is empty!" && exit 1
[ ! -f "a.MarkDuplicates.metrics.txt" ] && echo "MarkDuplicates output metrics file does not exist!" && exit 1
[ ! -s "a.MarkDuplicates.metrics.txt" ] && echo "MarkDuplicates output metrics file is empty!" && exit 1

echo "All tests succeeded!"
exit 0