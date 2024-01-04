#!/bin/bash

set -eo pipefail 

read_distribution.py \
    -i $par_input \
    -r $par_refgene \
> $par_output

# Version
read -r -d '' text <<- END_VERSIONS
"${meta_functionality_name}":
    rseqc: \$(read_distribution.py --version | sed -e "s/read_distribution.py //g")
END_VERSIONS

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
else
    echo "$text" > "$par_versions"
fi