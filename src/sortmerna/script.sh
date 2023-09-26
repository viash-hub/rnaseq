#!/bin/bash

set -eo pipefail

# IFS="," read -ra input <<< "$par_input"
# read_count="${#input[@]}"
read1=$(find $par_input/ -iname *_1_trimmed.fq.gz*)
read2=$(find $par_input/ -iname *_2_trimmed.fq.gz*)

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

    mv non_rRNA_reads.f*q.gz $par_output_1
    # mv rRNA_reads.log ${par_id}.sortmerna.log
else   
    sortmerna $refs \
    --reads "$read1" --reads "$read2" \
    --threads $meta_cpus --workdir . \
    --aligned rRNA_reads --fastx \
    --other non_rRNA_reads --paired_in --out2 
    # $args 
    
    mv non_rRNA_reads_fwd.f*q.gz ${par_output_1}
    mv non_rRNA_reads_rev.f*q.gz ${par_output_2}
    # mv rRNA_reads.log ${par_id}.sortmerna.log
fi