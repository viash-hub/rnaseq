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

python3 "$meta_resources_dir/tx2gene.py" \
    --quant_type $par_quant_type \
    --gtf $par_gtf \
    --quants $tmpdir \
    --id $par_gtf_group_features \
    --extra $par_gtf_extra_attributes \
    -o $par_tsv

# Version
# text="${meta_functionality_name}:
#     python: $(python3 --version | sed 's/Python //g')"

# if [ -e "$par_versions" ]; then
#     echo "$text" >> "$par_versions"
#     mv "$par_versions" "$par_updated_versions"
# else
#     echo "$text" > "$par_updated_versions"
# fi