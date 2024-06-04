#!/bin/bash

set -eo pipefail

filename="$(basename -- "$par_input")"

if [ ${filename##*.} == "gz" ]; then
    gunzip -c $par_input > $par_output
else
    cat $par_input > $par_output
fi

# # Version
# text="${meta_functionality_name}: 
#     gunzip: $(echo $(gunzip --version 2>&1) | sed -n 's/.* \([0-9]\+\.[0-9]\+\).*/\1/p')"
# if [ -e "$par_versions" ]; then
#     echo "$text" >> "$par_versions"
#     mv "$par_versions" "$par_updated_versions"
# else
#     echo "$text" > "$par_updated_versions"
# fi