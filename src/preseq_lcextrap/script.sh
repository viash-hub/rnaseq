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

text="${meta_functionality_name}:
    preseq: $(echo $(preseq 2>&1) | sed 's/^.*Version: //; s/Usage:.*\$//')"

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
    mv "$par_versions" "$par_updated_versions"
else
    echo "$text" > "$par_updated_versions"
fi