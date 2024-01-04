#!/bin/bash

set -eo pipefail

filename="$(basename -- "$par_input")"

if [ ${filename##*.} == "gz" ]; then
    gunzip -c $par_input > $par_output
else
    cat $par_input > $par_output
fi

# Version
read -r -d '' text <<- END_VERSIONS
"${meta_functionality_name}":
    gunzip: \$(echo \$(gunzip --version 2>&1) | sed 's/^.*(gzip) //; s/ Copyright.*\$//')
END_VERSIONS

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
else
    echo "$text" > "$par_versions"
fi