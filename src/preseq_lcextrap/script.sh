#!/bin/bash

set -eo pipefail

file=$(basename -- "$par_bam")
filename="${file%.*}"

samtools sort -o sorted_$filename.bam -n $par_bam
bedtools bamtobed -i sorted_$filename.bam > $filename.bed
bedtools sort -i $filename.bed > sorted_$filename.bed

preseq lc_extrap \
    sorted_$filename.bed \
    $par_extra_preseq_args \
    ${par_paired:+-pe} \
    -o $par_output

# Version
read -r -d '' text <<- END_VERSIONS
"${meta_functionality_name}":
    preseq: \$(echo \$(preseq 2>&1) | sed 's/^.*Version: //; s/Usage:.*\$//')
END_VERSIONS

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
else
    echo "$text" > "$par_versions"
fi