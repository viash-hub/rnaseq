#!/bin/bash

set -eo pipefail

filename="$(basename -- $par_fasta)"

python3 "$meta_resources_dir/filter_gtf_for_genes_in_genome.py" --gtf $par_gtf --fasta $par_fasta -o $par_filtered_gtf

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