#!/bin/bash

set -eo pipefail

IFS="," read -ra input <<< $par_input
n_fastq=${#input[@]}

required_args=("-p" "--probability" "-n" "--read-count")

for arg in "${required_args[@]}"; do
    if [[ "$par_extra_args" == *"$arg"* ]]; then
        echo "FQ/SUBSAMPLE requires --probability (-p) or --record-count (-n) specified in task.ext.args!"
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

# # Version
# text="${meta_functionality_name}:
#     fq: $(echo $(fq subsample --version | sed 's/fq-subsample //g'))"

# if [ -e "$par_versions" ]; then
#     echo "$text" >> "$par_versions"
#     mv "$par_versions" "$par_updated_versions"
# else
#     echo "$text" > "$par_updated_versions"
# fi