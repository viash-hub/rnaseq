#!/bin/bash

set -eo pipefail

eval umi_tools extract -I $par_input -S $par_output \
    --extract-method $par_umitools_extract_method \
    --bc-pattern $par_umitools_bc_pattern \
    --umi-separator $par_umitools_umi_separator

# if [ ! $par_input2 ]; then
#     # single end
#     eval umi_tools extract -I $par_input -S $par_output \
#     --extract-method $par_umitools_extract_method \
#     --bc-pattern $par_umitools_bc_pattern \
#     --umi-separator $par_umitools_umi_separator 
# else 
#     # paired end
#     eval umi_tools extract -I $$par_input --read2-in=$par_input2 \
#     -S "$par_id.umi_extract_1.fastq.gz" \
#     --read2-out "$par_id.umi_extract_2.fastq.gz" \
#     --extract-method $par_umitools_extract_method \
#     --bc-pattern $par_umitools_bc_pattern \
#     --bc-pattern2 $par_umitools_bc_pattern2 \
#     --umi-separator $umitools_umi_separator
# fi

# save_umi_intermeds

