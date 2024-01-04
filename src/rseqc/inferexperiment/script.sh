#!/bin/bash

set -eo pipefail 

infer_experiment.py \
    -i $par_input \
    -r $par_refgene \
    -s $par_sample_size \
    -q $par_map_qual \
> $par_output

# Version
read -r -d '' text <<- END_VERSIONS
"${meta_functionality_name}":
    rseqc: \$(infer_experiment.py --version | sed -e "s/infer_experiment.py //g")
END_VERSIONS

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
else
    echo "$text" > "$par_versions"
fi