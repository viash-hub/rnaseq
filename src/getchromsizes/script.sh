#!/bin/bash

set -eo pipefail
par_fasta="/home/nirmayi/data_intuitive/rnaseq.vsh/testData/test_output/ref.cat_additional.fasta_output"
par_fai="/home/nirmayi/data_intuitive/rnaseq.vsh/testData/test_output/fai"
par_sizes="/home/nirmayi/data_intuitive/rnaseq.vsh/testData/test_output/sizes"

samtools faidx $par_fasta -o $par_fai | cut -f 1,2 "$par_fasta.fai" > "$par_sizes"
cat "$par_fasta.fai" > $par_fai