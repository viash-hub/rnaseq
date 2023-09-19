#!/bin/bash

set -eo pipefail

IFS="," read -ra input <<< "$par_input"

count="${#input[@]}"

if [ "$par_paired" == "true" ]; then
    echo "Paired - $count"
    if [ "$count" -ne 2 ]; then
        echo "Paired end input requires two files"
        exit 1
    else
        trim_galore $par_extra_trimgalore_args --gzip "${input[*]}" -o $par_output
    fi
else
    echo "Not Paired - $count"
    if [ "$count" -ne 1 ]; then
        echo "Single end input requires one file"
        exit 1
    else
        trim_galore $par_extra_trimgalore_args --gzip $par_input -o $par_output
    fi
fi 