#!/bin/bash

# define input and output for script
input="salmon.tar.gz"
output="salmon"

# run executable and tests
echo "> Running $meta_functionality_name."

"$meta_executable" \
    --input "$meta_resources_dir/$input" \
    --output "$output"

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> Checking whether output can be found and has content"

[ ! -d "$output" ] && echo "$output file missing!" && exit 1
[ -z "$(ls -A $output)" ] && echo "$output is empty!" && exit 1

exit 0