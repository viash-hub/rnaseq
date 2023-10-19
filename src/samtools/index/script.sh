#!/bin/bash

set -eo pipefail

## VIASH START
meta_cpus=2
## VIASH END

mkdir -p $par_output
 filename="$(basename -- $par_input/*.bam)"
if $par_bam_csi_index; then
    samtools index -@ $meta_cpus --csi -o "$par_output/$filename.csi" $par_input/$filename
else
    samtools index -@ $meta_cpus --csi -o "$par_output/$filename.bai" $par_input/$filename
fi

cp $par_input/$filename $par_output/$filename
