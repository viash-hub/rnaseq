#!/bin/bash

set -eo pipefail

file=$(basename -- "$par_bam")
filename="${file%.*}"

samtools sort -o sorted_$filename.bam -n $par_bam
bedtools bamtobed -i sorted_$filename.bam > $filename.bed
bedtools sort -i $filename.bed > sorted_$filename.bed

ls -la

preseq lc_extrap \
    sorted_$filename.bed \
    $par_extra_preseq_args \
    ${par_paired:+-pe} \
    -B \
    -v \
    -o $par_output
     