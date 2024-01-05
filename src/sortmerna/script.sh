#!/bin/bash

set -eo pipefail

IFS="," read -ra input <<< "$par_input"

refs=()
while IFS="," read -r path 
do
    refs+=("--ref $path")
done < $par_ribo_database_manifest

if [ "$par_paired" == "false" ]; then
    sortmerna $refs \
    --reads "${input[0]}" \
    --threads $meta_cpus --workdir . \
    --aligned rRNA_reads --fastx --other non_rRNA_reads \
    # $args
    mv non_rRNA_reads.f*q.gz "$par_fastq_1"
else   
    sortmerna $refs \
    --reads ${input[0]} --reads ${input[1]} \
    --threads $meta_cpus --workdir . \
    --aligned rRNA_reads --fastx \
    --other non_rRNA_reads --paired_in --out2 
    # $args 
    mv non_rRNA_reads_fwd.f*q.gz $par_fastq_1
    mv non_rRNA_reads_rev.f*q.gz $par_fastq_2
fi

mv rRNA_reads.log $par_sortmerna_log

# Version
text="${meta_functionality_name}:
    sortmerna: $(echo $(sortmerna --version 2>&1) | sed 's/^.*SortMeRNA version //; s/ Build Date.*\$//')"

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
    mv "$par_versions" "$par_updated_versions"
else
    echo "$text" > "$par_updated_versions"
fi