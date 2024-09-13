#!/bin/bash

set -eo pipefail

IFS=";" read -ra input <<< $par_input
n_fastq=${#input[@]}

required_args=("-p" "--probability" "-n" "--read-count")
for arg in "${required_args[@]}"; do
    if [[ "$par_extra_args" == *"$arg"* ]]; then
        echo "FQ/SUBSAMPLE requires either --probability (-p) or --record-count (-n) to be specified with --extra_args"
        exit 1
    fi
done

if [ $n_fastq -eq 1 ]; then
    fq subsample $par_extra_args ${input[*]} --r1-dst $par_output_1
elif [ $n_fastq -eq 2 ]; then
    fq subsample $par_extra_args ${input[*]} --r1-dst $par_output_1 --r2-dst $par_output_2
else 
    echo "FQ/SUBSAMPLE only accepts 1 or 2 FASTQ files!"
    exit 1
fi
