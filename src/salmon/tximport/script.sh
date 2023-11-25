#!/bin/bash

set -eo pipefail

function clean_up {
    rm -rf "$salmon_tmpdir"
}
trap clean_up EXIT

salmon_tmpdir=$(mktemp -d "$meta_temp_dir/$meta_functionality_name-XXXXXXXX")

IFS="," read -ra salmon_results <<< $par_salmon_quant_results
for result in ${salmon_results[*]}
do 
    cp -r $result $salmon_tmpdir
done

Rscript "$meta_resources_dir/salmon_tximport.r" NULL $salmon_tmpdir salmon.merged $par_tx2gene_tsv
