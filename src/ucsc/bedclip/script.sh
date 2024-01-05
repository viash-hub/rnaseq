#!/bin/bash

set -eo pipefail

bedClip $par_input_bedgraph $par_sizes $par_output_bedgraph

# # Version
# text="${meta_functionality_name}:
#     ucsc: version?"

# if [ -e "$par_versions" ]; then
#     echo "$text" >> "$par_versions"
    # mv "$par_versions" "$par_updated_versions"
# else
#     echo "$text" > "$par_updated_versions"
# fi