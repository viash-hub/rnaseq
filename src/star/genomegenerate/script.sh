#!/bin/bash

set -eo pipefail

## VIASH START
meta_cpus=8
## VIASH END

filename="$(basename -- $par_fasta/*)"

samtools faidx "$par_fasta/$filename"
NUM_BASES=`gawk '{sum = sum + \$2}END{if ((log(sum)/log(2))/2 - 1 > 14) {printf "%.0f", 14} else {printf "%.0f", (log(sum)/log(2))/2 - 1}}' $par_fasta/$filename.fai`
mkdir $par_star_index
STAR --runMode genomeGenerate --genomeDir $par_star_index --genomeFastaFiles $par_fasta/$filename --sjdbGTFfile $par_gtf/* --runThreadN $meta_cpus  --genomeSAindexNbases $NUM_BASES 
