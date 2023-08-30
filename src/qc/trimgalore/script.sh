#!/bin/bash

set -eo pipefail

if [ -f $par_umi_extract_output_r2 ]; then
    # paired-end
    eval trim_galore $par_extra_trimgalore_args --paired --gzip $par_umi_extract_output $par_umi_extract_r2_output -o $par_output # --cores $par_cores
else
    # single-end
    eval trim_galore $par_extra_trimgalore_args --gzip $par_umi_extract_output -o $par_output # --cores $par_cores
fi