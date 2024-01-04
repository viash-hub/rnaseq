#!/bin/bash

set -eo pipefail

multiqc \
    -f -p \
    ${par_multiqc_title:+--title $par_multiqc_title} \
    --config $par_multiqc_custom_config \
    $par_input

# cat <<-END_VERSIONS > versions.yml
# "$key":
#     multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
# END_VERSIONS
