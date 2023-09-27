#!/bin/bash

set -eo pipefail

mkdir -p $par_output
read1=$(find $par_input/ -iname primary_1*)
read2=$(find $par_input/ -iname primary_2*)

refs=()
while IFS="," read -r path 
do
    refs+=("--ref $path")
done < $par_ribo_database_manifest

if [ "$par_paired" == "false" ]; then
    sortmerna $refs \
    --reads "$read1" \
    --threads $meta_cpus --workdir . \
    --aligned rRNA_reads --fastx --other non_rRNA_reads \
    # $args

    mv non_rRNA_reads.f*q.gz "$par_output/non_rRNA_reads.f*q.gz"
    # mv rRNA_reads.log ${par_id}.sortmerna.log
else   
    sortmerna $refs \
    --reads "$read1" --reads "$read2" \
    --threads $meta_cpus --workdir . \
    --aligned rRNA_reads --fastx \
    --other non_rRNA_reads --paired_in --out2 
    # $args 
    
    mv non_rRNA_reads_fwd.f*q.gz $par_output/non_rRNA_reads_1.fq.gz
    mv non_rRNA_reads_rev.f*q.gz $par_output/non_rRNA_reads_2.fq.gz
    # mv rRNA_reads.log ${par_id}.sortmerna.log
fi