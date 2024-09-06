#!/bin/bash

set -eo pipefail

cut -f 1,7 $par_biocounts | tail -n +3 | cat $par_biotypes_header - >> $par_featurecounts_multiqc

python3 "$meta_resources_dir/mqc_features_stat.py" \
    $par_featurecounts_multiqc \
    -s $par_id \
    -f rRNA \
    -o $par_featurecounts_rrna_multiqc
