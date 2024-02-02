#!/bin/bash

set -eo pipefail

mkdir -p $par_output_dir

qualimap rnaseq \
    --java-mem-size=$par_java_memory_size \
    --algorithm $par_algorithm \
    --num-pr-bases $par_pr_bases \
    --num-tr-bias $par_tr_bias \
    --sequencing-protocol $par_sequencing_protocol \
    -bam $par_input \
    -gtf $par_gtf \
    ${par_paired:+-pe} \
    ${par_sorted:+-s} \
    -outdir $par_output_dir \
    -outformat $par_output_format 

# Version

text="${meta_functionality_name}:
    qualimap: $(echo $(qualimap --help 2>&1) | grep -oP 'QualiMap v\.\K\d+\.\d+')"

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
    mv "$par_versions" "$par_updated_versions"
else
    echo "$text" > "$par_updated_versions"
fi
