#!/bin/bash

echo ">>> Teatinf $meta_functionality_name"

"$meta_executable" \
  --id test \
  --biocounts "$meta_resources_dir/test.featureCounts.txt" \
  --featurecounts_multiqc "test.biotype_counts_rrna_mqc.tsv"

echo ">>> Check whether output exists"
[ ! -f test.biotype_counts_rrna_mqc.tsv ] && echo "MultiQC input file for biotype QC does not exist!" && exit 1

echo "All tests succeeded!"
exit 0