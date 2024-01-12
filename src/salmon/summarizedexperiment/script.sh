#!/bin/bash

set -eo pipefail

mkdir -p $par_output

Rscript "$meta_resources_dir/salmon_summarizedexperiment.r" NULL $par_counts_gene $par_tpm_gene $par_tx2gene_tsv 

Rscript "$meta_resources_dir/salmon_summarizedexperiment.r" NULL $par_counts_gene_length_scaled $par_tpm_gene $par_tx2gene_tsv 

Rscript "$meta_resources_dir/salmon_summarizedexperiment.r" NULL $par_counts_gene_scaled $par_tpm_gene $par_tx2gene_tsv 

Rscript "$meta_resources_dir/salmon_summarizedexperiment.r" NULL $par_counts_transcript $par_tpm_transcript $par_tx2gene_tsv 

mv ${par_counts_gene%.*}.rds $par_output/
mv ${par_counts_gene_length_scaled%.*}.rds $par_output/
mv ${par_counts_gene_scaled%.*}.rds $par_output/
mv ${par_counts_transcript%.*}.rds $par_output/

# Version
summarizedExperiment_ver=$(Rscript -e "library(SummarizedExperiment); cat(as.character(packageVersion('SummarizedExperiment')))")
text="${meta_functionality_name}:
    r-base: $(echo $(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    bioconductor-summarizedexperiment: ${summarizedExperiment_ver}"

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
    mv "$par_versions" "$par_updated_versions"
else
    echo "$text" > "$par_updated_versions"
fi