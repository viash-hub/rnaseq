#!/bin/bash

set -eo pipefail

"$meta_resources_dir/salmon_tx2gene.py" --gtf $par_gtf --salmon $par_salmon_quant_results --id $par_gtf_group_features --extra $par_gtf_extra_attributes -o $par_tsv