#!/bin/bash

set -eo pipefail

prefix_forward="forward"
prefix_reverse="reverse"
if [ $par_strandedness == 'reverse' ]; then
    prefix_forward="reverse"
    prefix_reverse="forward"
fi

bedtools genomecov \
    -ibam $par_bam \
    -bg \
    -strand + \
    $par_extra_bedtools_args \
    | bedtools sort > $prefix_forward.bedGraph

bedtools genomecov \
    -ibam $par_bam \
    -bg \
    -strand - \
    $par_extra_bedtools_args \
    | bedtools sort > $prefix_reverse.bedGraph

mv $prefix_forward.bedGraph $par_bedgraph_forward
mv $prefix_reverse.bedGraph $par_bedgraph_reverse