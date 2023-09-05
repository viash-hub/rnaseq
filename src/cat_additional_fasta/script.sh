#!/bin/bash

set -eo pipefail

genome_name="${par_fasta##*/}"
genome_name="${genome_name%%.*}"

biotype_name=''
if [ biotype ]; then 
    biotype_name="-b $par_biotype"
fi

add_name="${par_additional_fasta##*/}"
add_name="${par_additional_fasta%%.*}"

name="${genome_name}_${add_name}"

# Use fasta2gtf.py to generate a GTF annotation file from the FASTA file
fasta2gtf.py -o "$add_name.gtf" $biotype_name $par_additional_fasta

cat $par_fasta $par_additional_fasta > $par_fasta_output
cat $par_gtf "$add_name.gtf" > $par_gtf_output