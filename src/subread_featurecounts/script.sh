#!/bin/bash

set -eo pipefail

if [ $par_strandedness == 'forward' ]; then
    strandedness=1
elif [ $par_strandedness == 'reverse' ]; then
    strandedness=2
else
    strandedness=0
fi

if $par_gencode; then 
    par_extra_featurecounts_args+=" -g gene_type"
else
    par_extra_featurecounts_args+=" -g $par_featurecounts_group_type"
fi

if $par_paired; then
    paired=' -p'
fi

par_extra_featurecounts_args+=" -t $par_featurecounts_feature_type"

featureCounts \
    $par_extra_featurecounts_args \
    $paired \
    -T $meta_cpus \
    -a $par_gtf \
    -s $strandedness \
    -o $par_counts \
    $par_bam

