#!/bin/bash

set -eo pipefail

mkdir -p $par_output

## Ensures --strip-components only applied when top level of tar contents is a directory
## If just files or multiple directories, place all in prefix
if [ ${par_input#*.} == "tar.gz" ]; then
    if [[ $(tar -taf $par_input | grep -o -P "^.*?\\/" | uniq | wc -l) -eq 1 ]]; then
        tar -C $par_output --strip-components 1 -xavf $par_input --no-same-owner
    else
        tar -C $par_output -xavf $par_input --no-same-owner
    fi
else
    cat $par_input > $par_output
fi

# Version
text="${meta_functionality_name}:
    untar: $(echo $(tar --version 2>&1) | grep -oP 'tar \(GNU tar\) \K\d+\.\d+')"

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
    mv "$par_versions" "$par_updated_versions"
else
    echo "$text" > "$par_updated_versions"
fi