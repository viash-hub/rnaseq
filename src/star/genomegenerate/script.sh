#!/bin/bash

set -eo pipefail

samtools faidx "$par_fasta"
NUM_BASES=`gawk '{sum = sum + \$2}END{if ((log(sum)/log(2))/2 - 1 > 14) {printf "%.0f", 14} else {printf "%.0f", (log(sum)/log(2))/2 - 1}}' $par_fasta.fai`
mkdir $par_star_index
STAR --runMode genomeGenerate --genomeDir $par_star_index --genomeFastaFiles $par_fasta --sjdbGTFfile $par_gtf --runThreadN $meta_cpus  --genomeSAindexNbases $NUM_BASES 

# Version
read -r -d '' text <<- END_VERSIONS
"${meta_functionality_name}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
END_VERSIONS

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
else
    echo "$text" > "$par_versions"
fi