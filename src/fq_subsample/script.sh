#!/bin/bash

set -eo pipefail

IFS="," read -ra input <<< $par_input
n_fastq=${#input[@]}

if [ ! $par_extra_args == *'-p'* ] && [ ! $par_extra_args == *'--probability'* ] && [ ! $par_extra_args == *'-n'* ] && [ ! $par_extra_args == *'--record-count'* ]; then
    echo "FQ/SUBSAMPLE requires --probability (-p) or --record-count (-n) specified in task.ext.args!"
    exit 1
fi

if [ $n_fastq -eq 1 ]; then
    fq subsample $par_extra_args ${input[*]} --r1-dst $par_output_1
elif [ $n_fastq -eq 2 ]; then
    fq subsample $par_extra_args ${input[*]} --r1-dst $par_output_1 --r2-dst $par_output_2
else 
    echo "FQ/SUBSAMPLE only accepts 1 or 2 FASTQ files!"
    exit 1
fi