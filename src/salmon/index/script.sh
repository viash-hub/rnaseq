#!/bin/bash

set -eo pipefail

## VIASH START
meta_cpus=8
# par_genome_fasta="/home/nirmayi/data_intuitive/rnaseq.vsh/testData/reference/genome.fasta"
# transcript_fasta="/home/nirmayi/data_intuitive/rnaseq.vsh/testData/reference/transcriptome.fasta"
## VIASH END

grep '^>' $par_genome_fasta | cut -d ' ' -f 1 > decoys.txt
gentrome="gentrome.fa"
if [ ${par_genome_fasta##*.} == "gz" ]; then
    grep '^>' <(gunzip -c $par_genome_fasta) | cut -d ' ' -f 1 > decoys.txt
    gentrome="gentrome.fa.gz"
fi

sed -i.bak -e 's/>//g' decoys.txt
cat $transcript_fasta $genome_fasta > $gentrome

salmon index --threads $meta_cpus -t $gentrome --skipQuant -i salmon