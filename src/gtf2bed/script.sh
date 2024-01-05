#!/bin/bash

set -eo pipefail 

perl "$meta_resources_dir/gtf2bed.pl" $par_gtf > $par_bed_output

# Version

text="${meta_functionality_name}:
    perl: $(echo $(perl --version 2>&1) | sed 's/.*v\\(.*\\)) built.*/\\1/')"

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
    mv "$par_versions" "$par_updated_versions"
else
    echo "$text" > "$par_updated_versions"
fi