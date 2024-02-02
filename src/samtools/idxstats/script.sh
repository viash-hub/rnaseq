#!/bin/bash

set -eo pipefail

samtools idxstats --threads ${meta_cpus:-1} $par_bam > $par_output

# Version
text="${meta_functionality_name}:
    samtools: $(echo $(samtools --version 2>&1) | grep -oP 'samtools \K\d+\.\d+')"

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
    mv "$par_versions" "$par_updated_versions"
else
    echo "$text" > "$par_updated_versions"
fi