#!/bin/bash

set -eo pipefail

function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT

tmpdir=$(mktemp -d "$meta_temp_dir/$meta_name-XXXXXXXX")

IFS=";" read -ra results <<< $par_quant_results
for result in ${results[*]}
do
    cp -r $result $tmpdir
done

python3 "$meta_resources_dir/tx2gene.py" \
    --quant_type $par_quant_type \
    --gtf $par_gtf \
    --quants $tmpdir \
    --id $par_gtf_group_features \
    --extra $par_gtf_extra_attributes \
    -o $par_tsv
