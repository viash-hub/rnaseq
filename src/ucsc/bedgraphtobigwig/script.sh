#!/bin/bash

set -eo pipefail

bedGraphToBigWig $par_bedgraph $par_sizes $par_bigwig

# # Version
# read -r -d '' text <<- END_VERSIONS
# "${meta_functionality_name}":
#     ucsc: version?
# END_VERSIONS

# if [ -e "$par_versions" ]; then
#     echo "$text" >> "$par_versions"
# else
#     echo "$text" > "$par_versions"
# fi