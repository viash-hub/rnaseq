#!/bin/bash

set -eo pipefail

python3 "$meta_resources_dir/prepare-for-rsem.py" \
    --stdin=$par_bam \
    --stdout=$par_output \
    --log=$par_log 
