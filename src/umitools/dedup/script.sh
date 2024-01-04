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

# Version
read -r -d '' text <<- END_VERSIONS
"${meta_functionality_name}":
    umitools: \$(umi_tools --version 2>&1 | sed 's/^.*UMI-tools version://; s/ *\$//')
END_VERSIONS

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
else
    echo "$text" > "$par_versions"
fi