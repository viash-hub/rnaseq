#!/bin/bash

set -eo pipefail

samtools flagstat --threads $meta_cpus $par_bam > $par_output

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