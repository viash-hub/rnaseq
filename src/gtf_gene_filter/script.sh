#!/bin/bash

set -eo pipefail

gtf_gene_filter/filter_gtf_for_genes_in_genome.py --gtf $par_gtf --fasta $par_fasta -o $par_filtered_gtf