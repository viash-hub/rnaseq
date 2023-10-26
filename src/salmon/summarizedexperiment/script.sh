#!/bin/bash

set -eo pipefail

mkdir -p $par_output

"$meta_resources_dir/salmon_summarizedexperiment.r" NULL $par_input/*gene_counts.tsv $par_input/*gene_tpm.tsv $par_tx2gene

"$meta_resources_dir/salmon_summarizedexperiment.r" NULL $par_input/*gene_counts_length_scaled.tsv $par_input/*gene_tpm.tsv $par_tx2gene

"$meta_resources_dir/salmon_summarizedexperiment.r" NULL $par_input/*gene_counts_scaled.tsv $par_input/*gene_tpm.tsv $par_tx2gene

"$meta_resources_dir/salmon_summarizedexperiment.r" NULL $par_input/*transcript_counts.tsv $par_input/*transcript_tpm.tsv $par_tx2gene

mv $par_input/*gene_counts.rds $par_output/
mv $par_input/*gene_counts_length_scaled.rds $par_output/
mv $par_input/*gene_counts_scaled.rds $par_output/
mv $par_input/*transcript_counts.rds $par_output/