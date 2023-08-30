#!/bin/bash

set -eo pipefail

# eval umi_tools extract -I $par_input -S $par_output \
#     --extract-method $par_umitools_extract_method \
#     --bc-pattern $par_umitools_bc_pattern \
#     --umi-separator $par_umitools_umi_separator

if [ -f $par_input_r2 ]; then
    # paired end
    eval umi_tools extract -I $par_input --read2-in=$par_input_r2 \
    -S $par_output --read2-out $par_r2_output \
    --extract-method $par_umitools_extract_method \
    --bc-pattern $par_umitools_bc_pattern \
    --bc-pattern2 $par_umitools_bc_pattern2 \
    --umi-separator $par_umitools_umi_separator
else 
    # single end
    eval umi_tools extract -I $par_input -S $par_output \
    --extract-method $par_umitools_extract_method \
    --bc-pattern $par_umitools_bc_pattern \
    --umi-separator $par_umitools_umi_separator 
fi

# save_umi_intermeds

