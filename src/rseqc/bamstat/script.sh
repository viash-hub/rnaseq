#!/bin/bash

set -eo pipefail 

bam_stat.py \
    --input $par_input \
    --mapq $par_map_qual \
> $par_output

# Version
read -r -d '' text <<- END_VERSIONS
"${meta_functionality_name}":
    rseqc: \$(bam_stat.py --version | sed -e "s/bam_stat.py //g")
END_VERSIONS

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
else
    echo "$text" > "$par_versions"
fi