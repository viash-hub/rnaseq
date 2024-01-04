#!/bin/bash

set -eo pipefail

if [ -f "$par_fasta" ]; then
    reference="--reference $par_fasta"
else
    reference=""
fi

samtools stats --threads $meta_cpus $reference $par_bam > $par_output

# Version
read -r -d '' text <<- END_VERSIONS
"${meta_functionality_name}":
    samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
END_VERSIONS

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
else
    echo "$text" > "$par_versions"
fi