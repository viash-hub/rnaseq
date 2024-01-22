#!/bin/bash

# define input and output for script
input="$meta_resources_dir/genes.gff.gz"
output="genes.gff"

# run executable and tests
echo "> Running $meta_functionality_name."

"$meta_executable" \
    --input "$input" \
    --output "$output"

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> Checking whether output can be found and has content"

[ ! -f "$output" ] && echo "$output file missing" && exit 1
[ ! -s "$output" ] && echo "$output file is empty" && exit 1

exit 0