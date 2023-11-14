#!/bin/bash

set -eo pipefail

multiqc_config="assets/multiqc_config.yml"

multiqc \
    -f \
    ${par_multiqc_title:+--title par_multiqc_title} \
    --config $par_multiqc_custom_config \
    .

# cat <<-END_VERSIONS > versions.yml
# "$key":
#     multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
# END_VERSIONS
