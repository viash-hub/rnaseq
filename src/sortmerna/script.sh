#!/bin/bash

set -eo pipefail

IFS="," read -ra input <<< "$par_input"
IFS=";" read -ra paths <<< "$par_ribo_database_manifest"
refs=""
for i in "${paths[@]}"
do
    refs+="-ref $i "
done

if [ "$par_paired" == "false" ]; then
    sortmerna \
        $refs \
        -reads ${input[0]} \
        --threads ${meta_cpus:-1} \
        --workdir . \
        --aligned rRNA_reads \
        --fastx \
        -num_alignments 1 \
        --other non_rRNA_reads
    mv non_rRNA_reads.f*q.gz "$par_fastq_1"
else   
    sortmerna \
        $refs \
        -reads ${input[0]} \
        --reads ${input[1]} \
        --threads ${meta_cpus:-1} \
        --workdir . \
        --aligned rRNA_reads \
        --fastx \
        --num_alignments 1 \
        --other non_rRNA_reads \
        --paired_in \
        --out2 
    mv non_rRNA_reads_fwd.f*q.gz $par_fastq_1
    mv non_rRNA_reads_rev.f*q.gz $par_fastq_2
fi

mv rRNA_reads.log $par_sortmerna_log
