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
read -r -d '' text <<- END_VERSIONS
"${meta_functionality_name}":
    qualimap: \$(echo \$(qualimap 2>&1) | sed 's/^.*QualiMap v.//; s/Built.*\$//')
END_VERSIONS

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
else
    echo "$text" > "$par_versions"
fi
