#!/bin/bash

set -eo pipefail

avail_mem=3072
if $par_only_build_index; then
    filename="$(basename -- $par_primary_ref/*)"
fi 

other_refs=()
while IFS="," read -r name path 
do
    other_refs+=("ref_$name=$path")
done < $par_bbsplit_fasta_list


if $par_only_build_index; then
    if [ -f "$par_primary_ref/$filename" ] && [ ${#other_refs[@]} -gt 0 ]; then
        bbsplit.sh -Xmx${avail_mem}M ref_primary=$par_primary_ref/$filename ${other_refs[@]} path=$par_bbsplit_index threads=$meta_cpus
    else
        echo "ERROR: Please specify as input a primary fasta file along with names and paths to non-primary fasta files."
    fi
else
    mkdir -p $par_filtered_output
    index_files=''
    if [ -d "$par_built_bbsplit_index" ]; then
    index_files="path=$par_built_bbsplit_index"
    elif [ -f "$par_primary_ref/$filename" ] && [ ${#other_refs[@]} -gt 0 ]; then
        index_files="ref_primary=$primary_ref/$filename ${other_refs[@]}"
    else
        echo "ERROR: Please either specify a BBSplit index as input or a primary fasta file along with names and paths to non-primary fasta files."
    fi
    if $par_paired; then
        read1=$(find $par_input/ -iname *_1*.fq.gz*)
        read2=$(find $par_input/ -iname *_2*.fq.gz*)
        fastq_in="in=$read1 in2=$read2"
        fastq_out="basename=${par_filtered_output}/%_#.fastq.gz"
    else
        read1=$(find $par_input/ -iname *trimmed.fq.gz*)
        fastq_in="in=$read1"
        fastq_out="basename=${par_filtered_output}/%.fastq.gz"
    fi
    bbsplit.sh -Xmx${avail_mem}M $index_files threads=$meta_cpus $fastq_in $fastq_out refstats=bbsplit_stats.txt
fi
