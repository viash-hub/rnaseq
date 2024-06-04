#!/bin/bash

set -eo pipefail 

read_distribution.py \
    -i $par_input \
    -r $par_refgene \
> $par_output

# Version
# text="${meta_functionality_name}:
#     rseqc: $(read_distribution.py --version | sed -e 's/read_distribution.py //g')"

# if [ -e "$par_versions" ]; then
#     echo "$text" >> "$par_versions"
#     mv "$par_versions" "$par_updated_versions"
# else
#     echo "$text" > "$par_updated_versions"
# fi