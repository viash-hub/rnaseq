#!/bin/bash

set -eo pipefail

multiqc \
    -f -p \
    ${par_multiqc_title:+--title $par_multiqc_title} \
    --config $par_multiqc_custom_config \
    $par_input

# Version
read -r -d '' text <<- END_VERSIONS
"${meta_functionality_name}":
    multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
END_VERSIONS

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
else
    echo "$text" > "$par_versions"
fi
