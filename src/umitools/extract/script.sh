#!/bin/bash

set -eo pipefail

function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT 

tmpdir=$(mktemp -d "$meta_temp_dir/$meta_functionality_name-XXXXXXXX")

IFS="," read -ra input <<< "$par_input"
IFS="," read -ra pattern <<< "$par_bc_pattern"

read_count="${#input[@]}"
pattern_count="${#pattern[@]}"

if [ "$par_paired" == "true" ]; then
    echo "Paired - Reads: $read_count bc_patterns: $pattern_count"
    if [ "$read_count" -ne 2 ] || [ "$pattern_count" -ne 2 ]; then
        echo "Paired end input requires two read files and two UMI patterns"
        exit 1
    else
        read1="$(basename -- ${input[0]})"
        read2="$(basename -- ${input[1]})"
        umi_tools extract \
            -I "${input[0]}" --read2-in="${input[1]}" \
            -S "$tmpdir/$read1" \
            --read2-out="$tmpdir/$read2" \
            --extract-method $par_umitools_extract_method \
            --bc-pattern "${pattern[0]}" \
            --bc-pattern2 "${pattern[1]}" \
            --umi-separator $par_umitools_umi_separator
        if [ $par_umi_discard_read == 1 ]; then
            # discard read 1
            cp $tmpdir/$read1 $par_fastq_1
        elif [ $par_umi_discard_read == 2 ]; then
            # discard read 2
            cp $tmpdir/$read2 $par_fastq_1
        else
            cp $tmpdir/$read1 $par_fastq_1
            cp $tmpdir/$read2 $par_fastq_2
        fi
    fi
else
    echo "Not Paired - $read_count"
    if [ "$read_count" -ne 1 ] || [ "$pattern_count" -ne 1 ]; then
        echo "Single end input requires one read file and one UMI pattern"
        exit 1
    else
        read1="$(basename -- ${input[0]})"
        umi_tools extract \
            -I "${input[0]}" -S "$tmpdir/$read1" \
            --extract-method $par_umitools_extract_method \
            --bc-pattern "${pattern[0]}" \
            --umi-separator $par_umitools_umi_separator 
        cp $tmpdir/$read1 $par_fastq_1
    fi
fi

# Version
read -r -d '' text <<- END_VERSIONS
"${meta_functionality_name}":
    umitools: \$(umi_tools --version 2>&1 | sed 's/^.*UMI-tools version://; s/ *\$//')
END_VERSIONS

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
else
    echo "$text" > "$par_versions"
fi