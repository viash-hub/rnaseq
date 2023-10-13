#!/bin/bash

## VIASH START
meta_resources_dir="..."
## VIASH END

set -eo pipefail

filename="$(basename -- $par_fasta/*)"
mkdir -p $par_filtered_gtf

"$meta_resources_dir/filter_gtf_for_genes_in_genome.py" --gtf $par_gtf/* --fasta $par_fasta/* -o $par_filtered_gtf/${filename%%.*}_genes.gtf