#!/bin/bash

set -eo pipefail

filename="$(basename -- $par_fasta)"

samtools faidx $par_fasta
cut -f 1,2 "$par_fasta.fai" > $par_sizes
mv "$par_fasta.fai" $par_fai

# Version
read -r -d '' text <<- END_VERSIONS
"${meta_functionality_name}":
    getchromsizes: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
END_VERSIONS

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
else
    echo "$text" > "$par_versions"
fi