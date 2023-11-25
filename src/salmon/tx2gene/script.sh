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

python3 "$meta_resources_dir/salmon_tx2gene.py" \
    --gtf $par_gtf \
    --salmon $salmon_tmpdir \
    --id $par_gtf_group_features \
    --extra $par_gtf_extra_attributes \
    -o $par_tsv