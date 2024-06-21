#!/bin/bash

set -eo pipefail

function clean_up {
    rm -rf "$salmon_tmpdir"
}
trap clean_up EXIT

tmpdir=$(mktemp -d "$meta_temp_dir/$meta_functionality_name-XXXXXXXX")

IFS="," read -ra results <<< $par_quant_results
for result in ${results[*]}
do 
    cp -r $result $tmpdir
done

Rscript "$meta_resources_dir/tximport.r" \
    NULL \
    $tmpdir \
    $par_quant_type.merged \
    $par_quant_type \
    $par_tx2gene_tsv
