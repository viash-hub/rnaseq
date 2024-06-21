#!/bin/bash

set -eo pipefail

IFS=";" read -ra read_1 <<< $par_read_1
IFS=";" read -ra read_2 <<< $par_read_2

filename=$(basename -- "${read_1[0]}")
if [ ${filename##*.} == "gz" ]; then
    command="zcat"
else
    command="cat"
fi

if [ ${#read_1[@]} -gt 0 ]; then
    $command ${read_1[*]} > $par_fastq_1
fi
if [ ${#read_2[@]} -gt 0 ]; then
    $command ${read_2[*]} > $par_fastq_2
fi
