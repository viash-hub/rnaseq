#!/bin/bash

set -eo pipefail

## VIASH START
meta_cpus=8
## VIASH END

filename="$(basename -- "$par_genome_fasta")"

grep '^>' $par_genome_fasta/$filename | cut -d ' ' -f 1 > decoys.txt
gentrome="gentrome.fa"
if [ ${filename##*.} == "gz" ]; then
    grep '^>' <(gunzip -c $par_genome_fasta/$filename) | cut -d ' ' -f 1 > decoys.txt
    gentrome="gentrome.fa.gz"
fi

sed -i.bak -e 's/>//g' decoys.txt
cat $par_transcriptome_fasta/* $par_genome_fasta/* > $gentrome

salmon index --threads $meta_cpus -t $gentrome -i $par_salmon_index