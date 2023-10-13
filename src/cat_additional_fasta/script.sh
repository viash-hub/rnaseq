#!/bin/bash

## VIASH START
meta_resources_dir="..."
## VIASH END

set -eo pipefail
genome_name="$(basename -- $par_fasta/*)"
add_name="$(basename -- $par_additional_fasta/*)"
mkdir -p $par_fasta_output
mkdir -p $par_gtf_output

# Use fasta2gtf.py to generate a GTF annotation file from the FASTA file
"$meta_resources_dir/fasta2gtf.py" \
  -o ${add_name%%.*}.gtf \
  ${par_biotype:+-b $par_biotype} \
  $biotype_name \
  $par_additional_fasta/*

cat $par_fasta/* $par_additional_fasta/* > $par_fasta_output/${genome_name%%.*}_${add_name%%.*}.fasta
cat $par_gtf/* ${add_name%%.*}.gtf > $par_gtf_output/${genome_name%%.*}_${add_name%%.*}.gtf