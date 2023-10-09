#!/bin/bash

set -eo pipefail

mkdir -p $par_output

args="--random-seed=100"

if $par_paired; then
    paired="--paired"
    args+=" --unpaired-reads=discard --chimeric-pairs=discard"
else
    paired=""
fi

if $par_get_output_stats; then
    stats="--output-stats $par_output/"
else
    stats=""
fi


PYTHONHASHSEED=0 umi_tools dedup -I $par_input/*.bam -S $par_output/deduplicated.bam $stats $paired $args