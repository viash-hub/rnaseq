#!/bin/bash

set -eo pipefail

function clean_up {
    rm -rf "$tmpdir"
}

tmpdir=$(mktemp -d "$meta_temp_dir/$meta_functionality_name-XXXXXXXX")

IFS="," read -ra input <<< "$par_input"
read_count="${#input[@]}"

if [ "$par_paired" == "true" ]; then
    echo "Paired - $read_count"
    if [ "$read_count" -ne 2 ]; then
        echo "Paired end input requires two read files and two UMI patterns"
        exit 1
    else
        trim_galore $par_extra_trimgalore_args --paired --gzip ${input[0]} ${input[1]} -o $tmpdir 
        read1=$(find $tmpdir/ -iname *_1*.fq.gz*)
        read2=$(find $tmpdir/ -iname *_2*.fq.gz*)
        cp $read1 $par_fastq_1
        cp $read2 $par_fastq_2
    fi
else
    echo "Not Paired - $read_count"
    if [ "$read_count" -ne 1 ]; then
        echo "Single end input requires one read file and one UMI pattern"
        exit 1
    else
        read1="$(basename -- ${input[0]})"
        trim_galore $par_extra_trimgalore_args --gzip ${input[0]} -o $tmpdir
        read1=$(find $tmpdir/ -iname *trimmed.fq.gz*)
        cp $read1 $par_fastq_1
    fi
fi

trap clean_up EXIT 
