#!/bin/bash

set -eo pipefail

if [ $par_strandedness == 'forward' ]; then
    strandedness=1
elif [ $par_strandedness == 'reverse' ]; then
    strandedness=2
else
    strandedness=0
fi

featureCounts \
    $par_extra_featurecounts_args \
    ${par_paired:+-p} \
    -T $meta_cpus \
    -a $par_gtf \
    -s $strandedness \
    -o $par_counts \
    $par_bam

mv $par_counts.summary $par_summary