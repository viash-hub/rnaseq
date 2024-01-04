#!/bin/bash

echo ">>> Testing $meta_functionality_name"

gunzip $meta_resources_dir/genes.gtf.gz

"$meta_executable" \
    --strandedness reverse \
    --bam $meta_resources_dir/test.bam \
    --annotation_gtf $meta_resources_dir/genes.gtf \
    --transcript_gtf transcripts.gtf \
    --coverage_gtf coverage.gtf \
    --abundance  abundance.txt \
    --ballgown ballgown 
    
echo ">> Checking if the correct files are present"
[ ! -d ballgown ] && echo "Directory 'ballgown' does not exist!" && exit 1
[ ! -f transcripts.gtf ] && echo "File 'transcripts.gtf' does not exist!" && exit 1
[ ! -f coverage.gtf ] && echo "File 'coverage.gtf' does not exist!" && exit 1
[ ! -f abundance.txt ] && echo "File 'abundance.txt' does not exist!" && exit 1

echo ">>> Test finished successfully"
exit 0
