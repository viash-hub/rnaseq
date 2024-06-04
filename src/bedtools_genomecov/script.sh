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
    $par_extra_bedtools_args | bedtools sort > $prefix_forward.bedGraph

bedtools genomecov \
    -ibam $par_bam \
    -bg \
    -strand - \
    $par_extra_bedtools_args | bedtools sort > $prefix_reverse.bedGraph

mv $prefix_forward.bedGraph $par_bedgraph_forward
mv $prefix_reverse.bedGraph $par_bedgraph_reverse

# Version
# text="${meta_functionality_name}:
#     bedtools: $(bedtools --version | sed -e "s/bedtools\ v//g")"

# if [ -e "$par_versions" ]; then
#     echo "$text" >> "$par_versions"
#     mv "$par_versions" "$par_updated_versions"
# else
#     echo "$text" > "$par_updated_versions"
# fi