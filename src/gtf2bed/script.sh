#!/bin/bash

set -eo pipefail 

perl "$meta_resources_dir/gtf2bed.pl" $par_gtf > $par_bed_output

# Version
read -r -d '' text <<- END_VERSIONS
"${meta_functionality_name}":
    perl: \$(echo \$(perl --version 2>&1) | sed 's/.*v\\(.*\\)) built.*/\\1/')
END_VERSIONS

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
else
    echo "$text" > "$par_versions"
fi