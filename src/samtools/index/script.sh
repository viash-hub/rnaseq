#!/bin/bash

set -eo pipefail

if $par_bam_csi_index; then
    samtools index -@ $meta_cpus --csi -o $par_output_csi $par_input
else
    samtools index -@ $meta_cpus -o $par_output_bai $par_input
fi

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