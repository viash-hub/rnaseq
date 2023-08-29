#!/bin/bash

set -eo pipefail

# single-end
eval trim_galore $par_extra_trimgalore_args --gzip $par_umi_extract_output -o $par_output # --cores $par_cores

# paired-end
# eval trim_galore $par_extra_trimgalore_args --paired --gzip $par_input $par_input2 # --cores $par_cores