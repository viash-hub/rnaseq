#!/bin/bash

set -eo pipefail

cut -f 1,7 $par_biocounts | tail -n +3 | cat $par_biotypes_header - >> $par_featurecounts_multiqc

python3 "$meta_resources_dir/mqc_features_stat.py" \
    $par_featurecounts_multiqc \
    -s $par_id \
    -f rRNA \
    -o $par_featurecounts_rrna_multiqc

# Version

text="${meta_functionality_name}:
    python: $(python3 --version | sed 's/Python //g')"

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
    mv "$par_versions" "$par_updated_versions"
else
    echo "$text" > "$par_updated_versions"
fi