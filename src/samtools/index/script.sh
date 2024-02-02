#!/bin/bash

set -eo pipefail

if $par_bam_csi_index; then
    samtools index -@ ${meta_cpus:-1} --csi -o $par_output_csi $par_input
else
    samtools index -@ ${meta_cpus:-1} -o $par_output_bai $par_input
fi

# Version
text="${meta_functionality_name}:
    samtools: $(echo $(samtools --version 2>&1) | grep -oP 'samtools \K\d+\.\d+')"

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
    mv "$par_versions" "$par_updated_versions"
else
    echo "$text" > "$par_updated_versions"
fi