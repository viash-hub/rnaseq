#!/bin/bash

set -eo pipefail

python3 "$meta_resources_dir/prepare-for-rsem.py" \
    --stdin=$par_bam \
    --stdout=$par_output \
    --log=$par_log 

# Version
read -r -d '' text <<- END_VERSIONS
"${meta_functionality_name}":
    python: \$(python --version | sed 's/Python //g')
END_VERSIONS

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
else
    echo "$text" > "$par_versions"
fi