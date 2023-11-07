#!/bin/bash

set -eo pipefail

qualimap \
    rnaseq \
    --java-mem-size=$par_java_memory_size \
    --algorithm $par_algorithm \
    --num-pr-bases $par_pr_bases \
    --num-tr-bias $par_tr_bias \
    --sequencing-protocol $par_sequencing_protocol \
    -bam $par_input \
    -gtf $par_gtf \
    ${par_paired:+-pe} \
    ${par_sorted:+-s} \
    -outdir $meta_temp_dir \
    -outformat $par_output_format 
