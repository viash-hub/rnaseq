#!/bin/bash

set -eo pipefail

multiqc \
    -f -p \
    ${par_multiqc_title:+--title $par_multiqc_title} \
    --config $par_multiqc_custom_config \
    $par_input

# Version

text="${meta_functionality_name}:
    multiqc: $( multiqc --version | sed -e 's/multiqc, version //g' )"

if [ $par_versions ] && [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
    mv "$par_versions" "$par_updated_versions"
else
    echo "$text" > "$par_updated_versions"
fi
