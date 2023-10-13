#!/bin/bash

set -eo pipefail

filename="$(basename -- $par_fasta/*)"
mkdir -p $par_sizes
mkdir -p $par_fai

samtools faidx $par_fasta/$filename
cut -f 1,2 "$par_fasta/$filename.fai" > "$par_sizes/$filename.sizes"
mv "$par_fasta/$filename.fai" $par_fai/$filename.fai