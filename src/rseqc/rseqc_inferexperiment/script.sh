#!/bin/bash

set -eo pipefail 

infer_experiment.py \
    -i $par_input \
    -r $par_refgene \
    -s $par_sample_size \
    -q $par_map_qual \
> $par_output

# Version

# text="${meta_functionality_name}:
#     rseqc: $(infer_experiment.py --version | sed -e 's/infer_experiment.py //g')"

# if [ -e "$par_versions" ]; then
#     echo "$text" >> "$par_versions"
#     mv "$par_versions" "$par_updated_versions"
# else
#     echo "$text" > "$par_updated_versions"
# fi