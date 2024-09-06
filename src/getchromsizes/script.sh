#!/bin/bash

set -eo pipefail

filename="$(basename -- $par_fasta)"

samtools faidx $par_fasta
cut -f 1,2 "$par_fasta.fai" > $par_sizes
mv "$par_fasta.fai" $par_fai
