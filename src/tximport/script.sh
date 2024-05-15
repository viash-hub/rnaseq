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

# Version
tximeta_ver=$(Rscript -e "library(tximeta); cat(as.character(packageVersion('tximeta')))")
summarizedExperiment_ver=$(Rscript -e "library(SummarizedExperiment); cat(as.character(packageVersion('SummarizedExperiment')))")
tximport_ver=$(Rscript -e "library(tximport); cat(as.character(packageVersion('tximport')))")
text="${meta_functionality_name}:
    r-base: $(echo $(R --version 2>&1) | grep -oP 'R version \K\d+\.\d+\.\d+')
    bioconductor-tximeta: ${tximeta_ver}
    bioconductor-SummarizedExperiment: ${summarizedExperiment_ver}
    bioconductor-tximport: ${tximport_ver}"

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
    mv "$par_versions" "$par_updated_versions"
else
    echo "$text" > "$par_updated_versions"
fi