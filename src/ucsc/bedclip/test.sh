#!/bin/bash

echo ">>> Testing $meta_functionality_name"

"$meta_executable" \
    --input_bedgraph $meta_resources_dir/test.bedgraph \
    --sizes $meta_resources_dir/chrom_sizes \
    --output_bedgraph test.clipped.bedgraph  
    
echo ">> Checking if the correct files are present"
[ ! -f test.clipped.bedgraph ] && echo "File 'test.clipped.bedgraph' does not exist!" && exit 1

echo ">>> Test finished successfully"
exit 0
