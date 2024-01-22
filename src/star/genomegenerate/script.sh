#!/bin/bash

set -eo pipefail

samtools faidx "$par_fasta"
NUM_BASES=`gawk '{sum = sum + \$2}END{if ((log(sum)/log(2))/2 - 1 > 14) {printf "%.0f", 14} else {printf "%.0f", (log(sum)/log(2))/2 - 1}}' $par_fasta.fai`

STAR \
    --runMode genomeGenerate \
    --genomeDir $par_star_index \
    --genomeFastaFiles $par_fasta \
    --sjdbGTFfile $par_gtf \
    --runThreadN ${meta_cpus:-1} \
    --genomeSAindexNbases $NUM_BASES 

# Version
text="${meta_functionality_name}:
    star: $(STAR --version | sed -e "s/STAR_//g")
    samtools: $(echo $(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    gawk: $(echo $(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')"

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
    mv "$par_versions" "$par_updated_versions"
else
    echo "$text" > "$par_updated_versions"
fi