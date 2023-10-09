#!/bin/bash

set -eo pipefail

## VIASH START
par_input="test_sorted.bam"
par_output="test_index"
par_bam_csi_index=false
meta_cpus=2
## VIASH END

mkdir -p $par_output
cp $par_input $par_output/$par_input.bam

if $par_bam_csi_index; then
    samtools index -@ $meta_cpus --csi -o "$par_output/$par_input.csi" $par_input
else
    samtools index -@ $meta_cpus --csi -o "$par_output/$par_input.bai" $par_input
fi

