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
read -r -d '' text <<- END_VERSIONS
"${meta_functionality_name}":
    cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
END_VERSIONS

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
else
    echo "$text" > "$par_versions"
fi