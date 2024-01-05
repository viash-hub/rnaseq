#!/bin/bash

set -eo pipefail

gffread $par_input -o $par_output

# Version

text="${meta_functionality_name}:
    gffread: $(gffread --version 2>&1)"

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
    mv "$par_versions" "$par_updated_versions"
else
    echo "$text" > "$par_updated_versions"
fi