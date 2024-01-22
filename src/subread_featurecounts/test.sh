#!/bin/bash

echo ">>> Testing $meta_functionality_name"

"$meta_executable" \
    --paired false \
    --strandedness forward \
    --bam $meta_resources_dir/wgEncodeCaltechRnaSeqGm12878R1x75dAlignsRep2V2.bam \
    --gtf $meta_resources_dir/genes.gtf \
    --extra_featurecounts_args "-B -C" \
    --featurecounts_group_type gene_id \
    --featurecounts_feature_type exon \
    --gencode false \
    --counts featureCounts.txt \
    --summary featureCounts.txt.summary
    
echo ">> Checking if the correct files are present"
[ ! -f "featureCounts.txt" ] && echo "File 'featureCounts.txt' does not exist!" && exit 1
[ ! -f "featureCounts.txt.summary" ] && echo "File 'featureCounts.txt.summary' does not exist!" && exit 1
[ ! -s "featureCounts.txt" ] && echo "File 'featureCounts.txt' is empty!" && exit 1
[ ! -s "featureCounts.txt.summary" ] && echo "File 'featureCounts.txt.summary' is empty!" && exit 1

echo ">>> Test finished successfully"
exit 0
