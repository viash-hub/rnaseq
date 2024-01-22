#!/bin/bash

set -eo pipefail

if [ -f "$par_fasta" ]; then
    reference="--reference $par_fasta"
else
    reference=""
fi

samtools stats --threads ${meta_cpus:-1} $reference $par_bam > $par_output

# Version
text="${meta_functionality_name}:
    samtools: $(echo $(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')"

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
    mv "$par_versions" "$par_updated_versions"
else
    echo "$text" > "$par_updated_versions"
fi