#!/bin/bash

echo ">>> Testing $meta_functionality_name"

gunzip $meta_resources_dir/genes.gtf.gz

"$meta_executable" \
    --input $meta_resources_dir/genome.fasta \
    --gtf $meta_resources_dir/genes.gtf \
    --star_index  STAR_index \

echo ">> Checking if the correct files are present"
[ ! -d STAR_index ] && echo "Directory 'STAR_index' does not exist!" && exit 1

echo ">>> Test finished successfully"
exit 0
