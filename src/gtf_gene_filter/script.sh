#!/bin/bash

set -eo pipefail

filename="$(basename -- $par_fasta)"

python3 "$meta_resources_dir/filter_gtf_for_genes_in_genome.py" --gtf $par_gtf --fasta $par_fasta -o $par_filtered_gtf