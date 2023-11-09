#!/bin/bash

set -eo pipefail

if $par_bam_csi_index; then
    samtools index -@ $meta_cpus --csi -o $par_output_csi $par_input
else
    samtools index -@ $meta_cpus -o $par_output_bai $par_input
fi
