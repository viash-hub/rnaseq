#!/bin/bash

set -eo pipefail

cut -f 1,7 $par_biocounts | tail -n +3 | cat $par_biotypes_header - >> biotype_counts_mqc.tsv

python3 "$meta_resources_dir/mqc_features_stat.py" \
    biotype_counts_mqc.tsv \
    -s $par_id \
    -f rRNA \
    -o $par_featurecounts_multiqc

# Version
read -r -d '' text <<- END_VERSIONS
"${meta_functionality_name}":
    python: \$(python --version | sed 's/Python //g')
END_VERSIONS

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
else
    echo "$text" > "$par_versions"
fi