#!/bin/bash

set -eo pipefail

IFS="," read -ra input <<< "$par_input"
read_count="${#input[@]}"

refs=()
while IFS="," read -r path 
do
    refs+=("--ref $path")
done < $par_ribo_database_manifest

if [ "$par_paired" == "false" ]; then
    echo "Not Paired - $count"
    if [ "$read_count" -ne 1 || "$pattern_count" -ne 1 ]; then
        echo "Single end input requires one read file and one UMI pattern"
        exit 1
    else
        sortmerna $refs \
        --reads "${par_input[0]}" \
        --threads $meta_cpus --workdir . \
        --aligned rRNA_reads --fastx --other non_rRNA_reads \
        # $args

        mv non_rRNA_reads.f*q.gz $par_output_1
        # mv rRNA_reads.log ${par_id}.sortmerna.log
    fi
else
    echo "Paired - $count"
    if [ "$read_count" -ne 2 || "$pattern_count" -ne 2 ]; then
        echo "Paired end input requires two read files and two UMI patterns"
        exit 1
    else
        sortmerna $refs \
        --reads "${par_input[0]}" --reads "${par_input[1]}" \
        --threads $meat_cpus --workdir . \
        --aligned rRNA_reads --fastx \
        --other non_rRNA_reads --paired_in --out2 
        # $args 
        
        mv non_rRNA_reads_fwd.f*q.gz ${par_output_1}
        mv non_rRNA_reads_rev.f*q.gz ${par_output_2}
        # mv rRNA_reads.log ${par_id}.sortmerna.log
    fi
fi