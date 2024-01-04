#!/bin/bash

set -eo pipefail

samtools sort -@ $meta_cpus -o $par_output $par_input

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