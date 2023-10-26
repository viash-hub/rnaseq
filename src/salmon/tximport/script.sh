#!/bin/bash

set -eo pipefail

# function clean_up {
#     rm -rf "$tmpdir"
# }

# tmpdir=$(mktemp -d "$meta_temp_dir/$meta_functionality_name-XXXXXXXX")
mkdir -p $par_output
"$meta_resources_dir/salmon_tximport.r" NULL $par_salmon_quant_results $par_output/salmon.merged $par_tx2gene_tsv

# trap clean_up EXIT 