#!/bin/bash

set -eo pipefail

file=$(basename -- "$par_input")
filename="${file%.*}"

if [ "${file##*.}" == "bam" ]; then 
    samtools sort -o sorted_$filename.bam -n $par_input
    bedtools bamtobed -i sorted_$filename.bam > $filename.bed
    bedtools sort -i $filename.bed > sorted_$filename.bed
elif [ "${file##*.}" == "bed" ]; then
    bedtools sort -i $par_input > sorted_$filename.bed
else 
    echo "Invalid input file format!"
    exit 1
fi

if $par_paired; then
    paired="-pe"
else
    paired=""
fi

preseq lc_extrap \
    sorted_$filename.bed \
    $paired \
    $par_extra_preseq_args \
    -o $par_output
