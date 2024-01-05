#!/bin/bash

set -eo pipefail

samtools sort -@ $meta_cpus -o $par_output $par_input

# Version
text="${meta_functionality_name}:
    samtools: $(echo $(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')"

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
    mv "$par_versions" "$par_updated_versions"
else
    echo "$text" > "$par_updated_versions"
fi