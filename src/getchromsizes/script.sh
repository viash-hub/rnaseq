#!/bin/bash

set -eo pipefail

samtools faidx $par_fasta 
cut -f 1,2 "$par_fasta.fai" > "$par_sizes"
cat "$par_fasta.fai" > "$par_fai"
