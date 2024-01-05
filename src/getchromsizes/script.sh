#!/bin/bash

set -eo pipefail

filename="$(basename -- $par_fasta)"

samtools faidx $par_fasta
cut -f 1,2 "$par_fasta.fai" > $par_sizes
mv "$par_fasta.fai" $par_fai

# Version

text="${meta_functionality_name}:
    getchromsizes: $(echo $(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')"

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
    mv "$par_versions" "$par_updated_versions"
else
    echo "$text" > "$par_updated_versions"
fi