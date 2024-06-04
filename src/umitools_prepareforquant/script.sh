#!/bin/bash

set -eo pipefail

python3 "$meta_resources_dir/prepare-for-rsem.py" \
    --stdin=$par_bam \
    --stdout=$par_output \
    --log=$par_log 

# # Version
# text="${meta_functionality_name}:
#     python: $(python3 --version | sed 's/Python //g')"

# if [ -e "$par_versions" ]; then
#     echo "$text" >> "$par_versions"
#     mv "$par_versions" "$par_updated_versions"
# else
#     echo "$text" > "$par_updated_versions"
# fi