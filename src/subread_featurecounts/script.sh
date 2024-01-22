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
    paired='-p'
fi

par_extra_featurecounts_args+=" -t $par_featurecounts_feature_type"
echo "featureCounts $par_extra_featurecounts_args $paired -T ${meta_cpus:-1} -a $par_gtf -s $strandedness -o $par_counts $par_bam"
featureCounts \
    $par_extra_featurecounts_args \
    $paired \
    -T ${meta_cpus:-1} \
    -a $par_gtf \
    -s $strandedness \
    -o $par_counts \
    $par_bam

# Version
text="${meta_functionality_name}:
    subread: $( echo $(featureCounts -v 2>&1) | sed -e 's/featureCounts v//g')"

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
    mv "$par_versions" "$par_updated_versions"
else
    echo "$text" > "$par_updated_versions"
fi