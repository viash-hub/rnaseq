#!/bin/bash

set -eo pipefail
# mkdir -p $par_output

IFS="," read -ra input <<< "$par_input"
IFS="," read -ra pattern <<< "$par_bc_pattern"

read_count="${#input[@]}"
pattern_count="${#pattern[@]}"

if [ "$par_paired" == "true" ]; then
    echo "Paired - $count"
    if [ "$read_count" -ne 2 || "$pattern_count" -ne 2 ]; then
        echo "Paired end input requires two read files and two UMI patterns"
        exit 1
    else
        umi_tools extract -I "${input[0]}" --read2-in="${input[1]}" \
        -S "$par_output_1" \
        --read2-out="$par_output_2" \
        --extract-method $par_umitools_extract_method \
        --bc-pattern "${pattern[0]}" \
        --bc-pattern2 "${pattern[1]}" \
        --umi-separator $par_umitools_umi_separator
    fi
else
    echo "Not Paired - $count"
    if [ "$read_count" -ne 1 || "$pattern_count" -ne 1 ]; then
        echo "Single end input requires one read file and one UMI pattern"
        exit 1
    else
        umi_tools extract -I "${input[0]}" -S "$par_output_1" \
        --extract-method $par_umitools_extract_method \
        --bc-pattern "${pattern[0]}" \
        --umi-separator $par_umitools_umi_separator 
    fi
fi