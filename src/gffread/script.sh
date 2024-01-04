#!/bin/bash

set -eo pipefail

gffread $par_input -o $par_output

# Version
read -r -d '' text <<- END_VERSIONS
"${meta_functionality_name}":
    gffread: \$(gffread --version 2>&1)
END_VERSIONS

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
else
    echo "$text" > "$par_versions"
fi