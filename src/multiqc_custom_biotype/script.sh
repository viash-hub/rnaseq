#!/bin/bash

set -eo pipefail

cut -f 1,7 $par_biocounts | tail -n +3 | cat $par_biotypes_header - >> biotype_counts_mqc.tsv

"$meta_resources_dir/mqc_features_stat.py" \
    biotype_counts_mqc.tsv \
    -s $par_id \
    -f rRNA \
    -o $par_featurecounts_multiqc