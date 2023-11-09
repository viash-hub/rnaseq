#!/bin/bash

set -eo pipefail

args="--random-seed=100"

if $par_paired; then
    paired="--paired"
    args+=" --unpaired-reads=discard --chimeric-pairs=discard"
else
    paired=""
fi

if $par_get_output_stats; then
    mkdir -p $par_output_stats
    stats="--output-stats $par_output_stats/"
else
    stats=""
fi

PYTHONHASHSEED=0 umi_tools dedup -I $par_bam -S $par_output_bam $stats $paired $args