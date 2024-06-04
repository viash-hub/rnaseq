#!/bin/bash

set -eo pipefail 

bam_stat.py \
    --input $par_input \
    --mapq $par_map_qual \
> $par_output

# Version

# text="${meta_functionality_name}:
#     rseqc: $(bam_stat.py --version | sed -e 's/bam_stat.py //g')"

# if [ -e "$par_versions" ]; then
#     echo "$text" >> "$par_versions"
#     mv "$par_versions" "$par_updated_versions"
# else
#     echo "$text" > "$par_updated_versions"
# fi