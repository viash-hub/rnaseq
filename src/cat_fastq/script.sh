#!/bin/bash

set -eo pipefail

IFS=";" read -ra read_1 <<< $par_read_1
IFS=";" read -ra read_2 <<< $par_read_2

# cat $read_1 > $par_fastq_1
# cat $read_2 > $par_fastq_2

if [ ${#read_1[@]} -gt 0 ]; then
    cat $read_1 > $par_fastq_1
fi
if [ ${#read_2[@]} -gt 0 ]; then
    cat $read_2 > $par_fastq_2
fi

# Version

text="${meta_functionality_name}:
    cat: $(echo $(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')"

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
    mv "$par_versions" "$par_updated_versions"
else
    echo "$text" > "$par_updated_versions"
fi