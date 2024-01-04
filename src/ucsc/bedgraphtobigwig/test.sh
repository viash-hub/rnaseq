#!/bin/bash

echo ">>> Testing $meta_functionality_name"

"$meta_executable" \
    --bedgraph $meta_resources_dir/test.bedgraph \
    --sizes $meta_resources_dir/chrom_sizes \
    --bigwig test.bigwig  
    
echo ">> Checking if the correct files are present"
[ ! -f test.bigwig ] && echo "File 'test.bigwig' does not exist!" && exit 1

echo ">>> Test finished successfully"
exit 0
