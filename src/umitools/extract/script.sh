#!/bin/bash

set -eo pipefail
mkdir -p $par_output

IFS="," read -ra input <<< "$par_input"
IFS="," read -ra pattern <<< "$par_bc_pattern"

read_count="${#input[@]}"
pattern_count="${#pattern[@]}"

if [ "$par_paired" == "true" ]; then
    echo "Paired - $read_count"
    if [ "$read_count" -ne 2 ] || [ "$pattern_count" -ne 2 ]; then
        echo "Paired end input requires two read files and two UMI patterns"
        exit 1
    else
        read1="$(basename -- ${input[0]})"
        read2="$(basename -- ${input[1]})"
        umi_tools extract -I "${input[0]}" --read2-in="${input[1]}" \
        -S "$par_output/$read1" \
        --read2-out="$par_output/$read2" \
        --extract-method $par_umitools_extract_method \
        --bc-pattern "${pattern[0]}" \
        --bc-pattern2 "${pattern[1]}" \
        --umi-separator $par_umitools_umi_separator
        cp $par_output/$read1 $par_fastq_1
        cp $par_output/$read2 $par_fastq_2
    fi
else
    echo "Not Paired - $read_count"
    if [ "$read_count" -ne 1 ] || [ "$pattern_count" -ne 1 ]; then
        echo "Single end input requires one read file and one UMI pattern"
        exit 1
    else
        read1="$(basename -- ${input[0]})"
        umi_tools extract -I "${input[0]}" -S "$par_output/$read1" \
        --extract-method $par_umitools_extract_method \
        --bc-pattern "${pattern[0]}" \
        --umi-separator $par_umitools_umi_separator 
        cp $par_output/$read1 $par_fastq_1
    fi
fi