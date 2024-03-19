#!/bin/bash

set -eo pipefail

kallisto \
    index \
    ${par_pseudo_aligner_kmer_size:+-k $par_pseudo_aligner_kmer_size} \
    -i $par_kallisto_index \
    $par_transcriptome_fasta

# Version
text="${meta_functionality_name}:
    salmon: $(echo $(kallisto 2>&1) | grep -oP 'kallisto \K\d+\.\d+\.\d+')"
if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
    mv "$par_versions" "$par_updated_versions"
else
    echo "$text" > "$par_updated_versions"
fi