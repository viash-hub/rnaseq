#!/bin/bash

set -eo pipefail

filename="$(basename -- $par_fasta)"

python3 "$meta_resources_dir/filter_gtf_for_genes_in_genome.py" --gtf $par_gtf --fasta $par_fasta -o $par_filtered_gtf

# Version

text="${meta_functionality_name}:
    python: $(python3 --version | sed 's/Python //g')"

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
    mv "$par_versions" "$par_updated_versions"
else
    echo "$text" > "$par_updated_versions"
fi