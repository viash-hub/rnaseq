#!/bin/bash

# define input and output for script
input="hisat2.tar.gz"
output="hisat2"

# create temporary directory
tmpdir=$(mktemp -d "$meta_temp_dir/$meta_functionality_name-XXXXXXXX")
function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT

# run executable and tests
echo "> Running $meta_functionality_name."

"$meta_executable" \
    --input "$meta_resources_dir/$input" \
    --output "$tmpdir/$output"

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> Checking whether output can be found and has content"

[[ ! -f "$tmpdir"/$output ]] && echo "$output file missing" && exit 1
[[ -z "$tmpdir"/$output ]] && echo "$output file is empty" && exit 1

exit 0